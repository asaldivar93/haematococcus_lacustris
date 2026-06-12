"""Translate genes from reference to strains."""

import importlib.resources
import json
import re

from pathlib import Path

import cobra
import plotly.express as px
import plotly.graph_objects as go
import plotly.io as pio
import polars as pl

from tqdm import tqdm

with importlib.resources.open_text("labutils", "plot_template.json") as file:
  style = json.load(file)

#style["layout"]["colorway"] = random.shuffle(style["layout"]["colorway"])
pio.templates["paper"] = go.layout.Template(
    data=style["data"],
    layout=style["layout"],
)
pio.templates.default = "simple_white+paper"

BOOLS = r"(?:and|or|not)"

def clean_bool_expr(s: str) -> str:
    # 1) normalize whitespace
    s = re.sub(r"\s+", " ", s).strip()

    # 2) remove leading/trailing operators (e.g. "and id2", "id3 or")
    s = re.sub(rf"^(?:{BOOLS})\b\s+", "", s, flags=re.I)
    s = re.sub(rf"\s+\b(?:{BOOLS})$", "", s, flags=re.I)

    # 3) fix consecutive operators: keep the *last* one ("or and" -> "and")
    # repeat until stable in case of longer chains ("or and or" -> "or")
    prev = None
    while prev != s:
        prev = s
        s = re.sub(rf"\b({BOOLS})\b(?:\s+\b({BOOLS})\b)+", r"\2", s, flags=re.I)

    # 4) clean up any extra spaces again
    s = re.sub(r"\s+", " ", s).strip()
    return s

def update_gene_ids(model, gene_map):
    print(f"\nUpdating genes in {model.name}")
    for rxn in tqdm(model.reactions):
        rule = rxn.gene_reaction_rule
        if not rule:
            continue
        gpr = rxn.gpr
        for gene in gpr.genes:
            old_id = gene
            new_id = gene_map[old_id] if gene in gene_map else ""
            rule = re.sub(rf"\b{re.escape(old_id)}\b", new_id, rule)
            rule = clean_bool_expr(rule)

        rxn.gene_reaction_rule = rule
        rxn.update_genes_from_gpr()

    print(f"Total genes before: {len(model.genes)}")
    genes_to_remove = [gene.id for gene in model.genes
                               if not bool(gene.reactions)]
    print(f"Removing {len(genes_to_remove)} unused genes")
    for gene_id in genes_to_remove:
        model.genes.remove(gene_id)
    print(f"Total genes after update: {len(model.genes)}")

def get_homologs(rbh_lazy, bh_lazy, cutoff_expr):

    query_rbh = (
        rbh_lazy
        .filter(
            cutoff_expr,
            pl.col("query").is_in(genes_in_model)
        )
        .sort("query", "bits", descending=True)
        .unique(subset="query", keep="first")
        .with_columns(method=pl.lit("rbh"))
    )

    query_bh = (
        bh_lazy
        .filter(
            cutoff_expr,
            pl.col("query").is_in(genes_in_model)
        )
        .sort("query", "bits", descending=True)
        .unique(subset="query", keep="first")
        .with_columns(method=pl.lit("bh"))
        .select(pl.col("*").exclude("alln"))
    )
    concat = (
        pl.concat([query_rbh, query_bh])
        .sort("query", "bits", descending=True)
        .unique(subset="query", keep="first")
    )
    homologs = concat.collect()
    missing_genes = set(genes_in_model) - set(homologs["query"].to_list())

    return homologs, missing_genes

ID_RE = re.compile(r"[A-Za-z]+\d+\.\d+")

def eval_rule(rule: str, present: set[str]) -> bool:
    if rule == "NONE":
        return True
    expr = ID_RE.sub(lambda m: str(m.group(0) in present), rule)
    # assumes rule uses "and"/"or"/parentheses already
    return bool(eval(expr, {"__builtins__": {}}, {}))

if __name__=="__main__":
    hits_folder = Path("data/1_interim/genomes/")
    model_path = Path("models/draft/v0.0.3/hlacustris.xml")
    out_dir = Path("models/draft/v0.0.4")
    out_dir.mkdir(parents=True, exist_ok=True)
    model = cobra.io.read_sbml_model(model_path)
    genes_in_model = [gene.id for gene in model.genes]
    gpr_by_rxn = [{"rxn_id": rxn.id, "gpr": rxn.gene_reaction_rule}
                  for rxn in model.reactions
                  if rxn.gene_reaction_rule]
    gprs_df = pl.DataFrame(gpr_by_rxn)
    strains = ["nies", "haep"]
    rbh_dict = {}
    bh_dict = {}
    homologs_dict = {}
    missing_dict = {}
    models_dict = {}

    for s in strains:
        rbh_path = hits_folder / f"rbh_red_{s}.csv"
        bh_path = hits_folder / f"dict_{s}_to_red.csv"

        rbh_lazy = pl.scan_csv(
            rbh_path,
            separator = "\t",
            has_header = False,
            new_columns = ["query","target","pident","evalue","bits","qcov","tcov","qlen","tlen"]
        )
        bh_lazy = pl.scan_csv(
            bh_path,
            separator = "\t",
            has_header = False,
            new_columns = ["query", "target", "pident", "evalue", "bits","qcov","tcov","qlen","tlen", "alln"]
        )
        cutoff_expr = (
            pl.col("pident")>=0,
            pl.col("tcov")>=0,
            pl.col("qcov")>=0,
            pl.col("evalue")<=1e-5,
        )
        homologs, missing = get_homologs(rbh_lazy, bh_lazy, cutoff_expr)


        rbh_dict[s] = rbh_lazy
        bh_dict[s] = bh_lazy
        homologs_dict[s] = homologs
        missing_dict[s] = missing
        query = homologs_dict[s].filter(pl.col("method")=="rbh")
        gene_map = {row["query"]: row["target"]
                    for row in query.iter_rows(named=True)}
        genes_in_map = set(query["query"])
        new_gprs_df = gprs_df.with_columns(
            pl.col("gpr")
            .map_elements(
                lambda r: eval_rule(r, genes_in_map),
                return_dtype=pl.Boolean
            )
            .alias(s)
        )
        rxns_to_remove = new_gprs_df.filter(~pl.col(s))["rxn_id"].to_list()

        models_dict[s] = model.copy()
        models_dict[s].remove_reactions(rxns_to_remove)
        mets_to_remove = [met for met in models_dict[s].metabolites
                                 if not bool(met.reactions)]
        models_dict[s].remove_metabolites(mets_to_remove)
        models_dict[s].id = f"iHlct_{s}"
        models_dict[s].name = f"H. lacustris strain {s}"
        update_gene_ids(models_dict[s], gene_map)
        cobra.io.write_sbml_model(models_dict[s], out_dir / f"{s}.xml")

    ## Plot histograms of all proteins
    q = (
        bh_dict["nies"]
        .filter(cutoff_expr)
    )
    q = q.collect()

    fig = go.Figure()
    fig.add_trace(
        go.Violin(
            x=q["qlen"], name="red"
        )
    )
    fig.add_trace(
        go.Violin(
            x=q["tlen"], name="nies"
        )
    )
    q = (
        bh_dict["haep"]
        .filter(cutoff_expr)
    )
    q = q.collect()
    fig.add_trace(
        go.Violin(
            x=q["tlen"], name="haep"
        )
    )
    fig.update_layout({
        "xaxis": {"title": {"text": "protein len"}}
    })
    fig.write_image("figures/prot_len_hist.png")

    ## Plot histograms for proteins in model
    fig = go.Figure()
    fig.add_trace(
        go.Violin(
            x=homologs_dict["nies"]["qcov"], name="nies"
        )
    )
    fig.add_trace(
        go.Violin(
            x=homologs_dict["haep"]["qcov"], name="haep"
        )
    )
    fig.update_layout({
        "xaxis": {"title": {"text": "query cov"}}
    })
    fig.write_image("figures/prot_cov_hist.png")

    model.genes[2]

    model.reactions.HMR_2215.gpr.remove_gene()
