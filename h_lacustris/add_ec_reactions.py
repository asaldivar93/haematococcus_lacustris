"""Script to gather MetaNetX reactions from protein EC annotations."""

import argparse
import re
from pathlib import Path

import cobra
import polars as pl
from Bio import SeqIO
from cobra import Metabolite, Reaction
from xlsxwriter import Workbook

from h_lacustris.databases import (
    BIGG_METS_DB,
    BIGG_RXNS_DB,
    METANETX_METS_DB,
    METANETX_RXNS_DB,
    MTNTX_COMPARTMENTS,
)

LOCATION_REGEX = "(" + "|".join(map(re.escape, MTNTX_COMPARTMENTS)) + ")"


def query_ec_hits(file_path: Path) -> pl.LazyFrame:
    """Parse EC annotation files.

    - file_path:
        File with columns query, target, ec-code, location,
    evalue, pident, qcov, tcov
    """
    hits = pl.scan_csv(file_path)
    return (
        hits
        # Keep only significant hits with EC annotation
        .filter(
            pl.col("ec-code").is_not_null(),
            pl.col("evalue") <= eval_cut,
            pl.col("pident") >= pident_cut,
            pl.col("qcov") >= qcov_cut,
            pl.col("tcov") >= tcov_cut,
            ~pl.col("query").is_in(genes_in_model),
        )
        .with_columns(
            # Make strings into lists
            pl.col("ec-code").str.split(";"),
            pl.col("location").str.to_lowercase().str.split(";"),
            # Score as the weighted average of pident, qcov, tcov
            score=0.4 * (pl.col("pident") / 100)
            + 0.3 * pl.col("qcov")
            + 0.3 * pl.col("tcov"),
        )
        # replace location information with compartments keyword
        .with_columns(
            pl.col("location").list.eval(
                pl.element().str.extract(LOCATION_REGEX, group_index=1),
            ),
        )
        # Make it one line per query per ec-code and per location
        .explode("ec-code")
        .explode("location")
        # For annotations that are identical in ec-code and location,
        # get the mean score for all hits and gather information on the hits
        .group_by("query", "ec-code", "location")
        .agg(
            [
                pl.col("score").mean(),
                pl.struct(["target", "score"]),
            ],
        )
        .sort("query", "score", descending=True)
    )
    return hits_lazy


def match_ec_to_rxns(hits_lazy: pl.LazyFrame) -> pl.LazyFrame:
    """Get reactions from ec-code."""
    rxn_db = METANETX_RXNS_DB
    # Get reactions with a matching ec-code in swissprot hits
    matches = (
        rxn_db
        # Keep only reactions with ec-code annotation
        .filter(pl.col("ec-code").is_not_null())
        # Make it one ec-cod per row
        .explode("ec-code")
        # Join to hits on ec-code
        .join(hits_lazy, on="ec-code", how="inner")
    )

    return (
        matches
        # Bring together all genes that have a match
        # for the same reaction in the same compartment
        .with_columns(
            pl.col("location").is_null().alias("loc_null"),
            pl.col("location").fill_null("cytoplasm"),
        )
        .group_by("id", "location")
        # Get info on the target and quality of annotation
        .agg(
            [
                pl.col("score").mean(),
                pl.col("query").unique().alias("gpr"),
                pl.struct(["query", "target"]),
            ],
        )
        # Build the gpr
        .with_columns(
            pl.col("gpr").list.join(separator=" or "),
        )
        # Add reaction info
        .join(METANETX_RXNS_DB, on="id", how="inner")
        .join(BIGG_RXNS_DB, on="bigg.reaction", how="left")
        .filter(
            # Keep only reactions with bigg_id
            pl.col("bigg.reaction").is_not_null(),
            # Keep reactions with compartment
            pl.col("location").is_not_null(),
            # Remove Recon3 and mus musculus reactions
            ~pl.col("reaction_string").str.contains(
                r"(C\d{5}|CE\d{4}|HC\d{5}|M\d{5})",
            ),
        )
    )


def get_metabolites(rxns_lazy: pl.LazyFrame) -> pl.LazyFrame:
    """Get metabolites and annotation from a reactions database."""
    return (
        rxns_lazy
        # Get metabolite ids from reactions string
        .select(
            pl.col("reaction_string")
            # remove stoichiometric coefficients (e.g. "2.0 ")
            .str.replace_all(r"\b\d+(\.\d+)?\s+", "")
            # replace arrows with "+"
            .str.replace_all(r"\s*(<->|->|<-)\s*", " + ")
            # split into metabolites
            .str.split(" + ")
            .alias("met_id"),
        )
        .explode("met_id")
        # trim whitespace
        .with_columns(
            pl.col("met_id").str.strip_chars(),
        )
        # remove empty strings
        .filter(pl.col("met_id") != "")
        .unique()
        # Get universal Bigg id
        .with_columns(
            pl.col("met_id").str.head(-2).alias("universal_bigg_id"),
        )
        # Get Bigg annotation information
        .join(BIGG_METS_DB, on="universal_bigg_id", how="inner")
        # Get MetaNetX annotation information
        .sort("met_id")
        .unique(subset="met_id")
        .join(
            METANETX_METS_DB,
            left_on="metanetx.chemical",
            right_on="id",
            how="left",
        )
    )

def write_excel(model, rxns_df, file_name):
    new_rxns = rxns_df["bigg.reaction"].to_list()
    new_genes = rxns_df.with_columns(pl.col("gpr").str.split(" ")).explode("gpr")["gpr"].unique().to_list()
    reactions = []
    for rxn in model.reactions:
        is_new=False
        if rxn.id in new_rxns:
            is_new=True
        this_rxn = {
            "id": rxn.id,
            "name": rxn.name,
            "reaction": rxn.build_reaction_string(),
            "gpr": rxn.gpr.to_string(),
            "subsystem": rxn.subsystem,
            "added_from_sw": is_new,
        }
        this_rxn.update(rxn.annotation)
        reactions.append(this_rxn)

    metabolites = []
    for met in model.metabolites:
        this_met = {
            "id": met.id,
            "name": met.name,
            "formula": met.formula,
            "charge": met.charge,
        }
        this_met.update(met.annotation)
        metabolites.append(this_met)

    genes = []
    for gen in model.genes:
        is_new=False
        if gen.id in new_genes:
            is_new=True
        this_gene = {
            "id": gen.id,
            "added_from_sw": is_new,
        }
        genes.append(this_gene)

    with Workbook(file_name) as wb:
        pl.DataFrame(reactions).write_excel(
            workbook=wb,
            worksheet="reactions",
        )
        pl.DataFrame(metabolites).write_excel(
            workbook=wb,
            worksheet="metabolites",
        )
        pl.DataFrame(genes).write_excel(
            workbook=wb,
            worksheet="genes",
        )


def parser():
    """Comand line parser."""
    parser = argparse.ArgumentParser(
        description="Add reactions based on EC annotation.",
    )
    parser.add_argument(
        "inputs_list",
        type=str,
        help="Tab separated file of strains annotation_paths and model_paths.",
    )
    parser.add_argument(
        "out_dir",
        type=str,
        help="Output_dir",
    )
    parser.add_argument(
        "--eval",
        default="1e-5",
        type=float,
        help="Min e-value of the annotation.",
    )
    parser.add_argument(
        "--pident",
        default="40",
        type=float,
        help="Min percentage identity of the annotation.",
    )
    parser.add_argument(
        "--qcov",
        default="0.6",
        type=float,
        help="Min query coverage of the annotation.",
    )
    parser.add_argument(
        "--tcov",
        default="0.6",
        type=float,
        help="Min target coverage of the annotation.",
    )

    return parser.parse_args()


if __name__ == "__main__":
    # Inputs
    args = parser()
    inputs_path: str = args.inputs_list
    eval_cut: float = args.eval
    pident_cut: float = args.pident
    qcov_cut: float = args.qcov
    tcov_cut: float = args.tcov
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    eval_cut: float = 1e-5
    pident_cut: float = 40
    qcov_cut: float = 0.6
    tcov_cut: float = 0.6
    out_dir = Path("test/")
    out_dir.mkdir(parents=True, exist_ok=True)

    inputs_path="data/add_ec.list"
    worklist = pl.read_csv(inputs_path)
    for row in worklist.iter_rows(named=True):
        strain_dir = out_dir / f"{row["id"]}"
        strain_dir.mkdir(exist_ok=True, parents=True)
        strain = row["id"]
        model_path = Path(row["model"])
        hits_path = Path(row["annotation"])
        faa_path = Path(row["genomes"])

        # Load base model
        print(f"Loading model for strain {strain} in {model_path}")
        model = cobra.io.read_sbml_model(model_path)
        mets_in_model = [met.id for met in model.metabolites]
        biggs_in_model = [rxn.id for rxn in model.reactions]
        mtntxs_in_model = [
            rxn.annotation.get("metanetx.reaction", "")
            for rxn in model.reactions
        ]
        genes_in_model = [g.id for g in model.genes]

        # Load annotation
        hits_lazy = query_ec_hits(hits_path)
        result_lazy = match_ec_to_rxns(hits_lazy)

        # Get reactions not in the model already
        new_rxns = (
            result_lazy
            # Get the compartment of the reaction
            .with_columns(
                pl.col("reaction_string")
                .str.split(" + ")
                .list.first()
                .str.tail(1)
                .str.replace_many(model.compartments)
                .alias("compartment"),
            ).filter(
                # Keep metabolic reactions only
                pl.col("transport").is_null(),
                # Reactions not in the model already
                ~pl.col("bigg.reaction").is_in(biggs_in_model),
                ~pl.col("id").is_in(mtntxs_in_model),
                # Reactions which compartment match the hit location annotation
                pl.col("location") == pl.col("compartment"),
            )
        )

        metabolites = get_metabolites(new_rxns)
        new_mets = metabolites.filter(~pl.col("met_id").is_in(mets_in_model))

        new_rxns = new_rxns.collect()
        new_mets = new_mets.collect()
        print(f"Adding {new_rxns.shape} new reactions")
        print(f"Adding {new_mets.shape} new metabolites")

        # Add new metabolites to the model
        mets_to_add = []
        for row in new_mets.iter_rows(named=True):
            # Create new metabolite
            new_met = Metabolite(
                id=row["met_id"],
                name=row["name"],
                formula=row["formula"],
                charge=row["charge"],
                compartment=row["met_id"][-1],
            )
            # Get the annotation info
            # TO DO: replace this to match the required databases names
            annotation = dict(list(row.items())[-14:])
            # Parse the annotation info
            if all(val is None for val in annotation.values()):
                annotation = {}
            else:
                annotation = {
                    db: value
                    for db, value in annotation.items()
                    if value is not None
                }
            new_met.annotation = annotation
            # Append to metabolites list
            mets_to_add.append(new_met)
        # Add new mets to model
        model.add_metabolites(mets_to_add)

        select_columns = [
            "reaction_string",
            "name",
            "gpr",
            "bigg.reaction",
            "ec-code",
            "seed.reaction",
            "metacyc.reaction",
            "sabiork.reaction",
            "kegg.reaction",
        ]
        rxns_to_add = []
        for row in new_rxns.select(select_columns).iter_rows(named=True):
            # Create a new Reaction
            new_rxn = Reaction(
                id=row["bigg.reaction"],
                name=row["name"] if row["name"] is not None else "",
            )
            new_rxn.gene_reaction_rule = row["gpr"]
            # Get the annotation info
            # TO DO: replace this to match the required databases names
            annotation = dict(list(row.items())[-6:])
            if all([val is None for val in annotation.values()]):
                annotation = {}
            else:
                annotation = {
                    db: value
                    for db, value in annotation.items()
                    if value is not None
                }
            new_rxn.annotation = annotation
            # Append to rxns list
            rxns_to_add.append(new_rxn)
        # Add reactions from the model
        model.add_reactions(rxns_to_add)
        # Add stoichiometry to new reactions
        for row in new_rxns.iter_rows(named=True):
            rxn = model.reactions.get_by_id(row["bigg.reaction"])
            rxn.build_reaction_from_string(row["reaction_string"])

        # Save the model
        model.repair()
        model.id = strain
        model.name = strain
        print(f"Saving to {strain_dir}")
        cobra.io.write_sbml_model(model, strain_dir / f"{strain}.xml")
        file_name = strain_dir / f"{strain}.xlsx"
        write_excel(model, new_rxns, file_name)
        fasta_name = strain_dir / f"{strain}.faa"
        genes_in_model = [g.id for g in model.genes]
        with Path.open(fasta_name, "w") as out_handle:
            for rec in SeqIO.parse(faa_path, "fasta"):
                # rec.id is already "first word" of the header for FASTA in BioPython
                if rec.id in genes_in_model:
                    SeqIO.write(rec, out_handle, "fasta")
