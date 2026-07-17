from pathlib import Path

import cobra
import polars as pl
from labutils.cobra.io import write_excel

from h_lacustris.databases import (
    METANETX_METS_DB,
    METANETX_RXNS_DB,
)

databases = ["ec-code", "bigg.reaction", "seed.reaction", "metacyc.reaction",
    "kegg.reaction", "reactome", "rh"]
mtntx_db = (
    METANETX_RXNS_DB
    .collect()
    .select(["id", *databases])
    .rename({"rh": "rhea"})
    .with_columns(
        pl.col("ec-code").list.join("|")
    )
)

def update_metabolites(metabolites_df: pl.DataFrame) -> list[cobra.Metabolite]:
    """Create list of metabolite instances to add to a model."""
    mets_to_add = []
    for met in metabolites_df.iter_rows(named=True):
        # Get met info from the database
        new_met = cobra.Metabolite(
            id=met["met_id"],
            formula=met["formula"],
            name=met["name"],
            charge=met["charge"],
            compartment=met["met_id"][-1],
        )

        # Only the annotations should be left in the row dict
        db_xref = ["metanetx.chemical", "metacyc.compound", "inchikey"]
        annotation = {}
        for db in db_xref:
            if met[db] is not None:
                annotation[db] = met[db]

        new_met.annotation = annotation

        mets_to_add.append(new_met)

    return mets_to_add

def update_reactions(reactions_df: pl.DataFrame) -> list[cobra.Reaction]:
    """Create a list of reaction instances to add to the model."""
    rxns_to_add = []
    for row in reactions_df.iter_rows(named=True):
        # Get reaction info from the dataframe
        new_rxn = cobra.Reaction(
            id=row["rxn_id"],
            name=row["name"],
            subsystem=row["subsystem"],
            lower_bound=row["lower_bound"],
            upper_bound=row["upper_bound"],
        )
        if row["gpr"]:
            new_rxn.gene_reaction_rule = row["gpr"]

        # Only the annotations should be left in the row dict
        db_xref = ["metanetx.reaction", "ec-code", "metacyc.reaction"]
        annotation = {}
        for db in db_xref:
            if row[db] is not None:
                annotation[db] = row[db]
        new_rxn.annotation = annotation
        rxns_to_add.append(new_rxn)
    return rxns_to_add

def update_boundary_rxns(model, boundary):
    for row in boundary.iter_rows(named=True):
        met = model.metabolites.get_by_id(row["met_id"])
        met.compartment = met.id[-1]
        model.add_boundary(met, type=row["type"])
        if row["type"]=="exchange":
            rxn = model.reactions.get_by_id("EX_" + met.id)
            rxn.lower_bound=0
            rxn.upper_bound=0

def update_dels(model, deletions):
    mets_to_remove = []
    rxns_to_remove = []
    for row in deletions.iter_rows(named=True):
        match row["type"]:
            case "metabolite":
                met = model.metabolites.get_by_id(row["id"])
                mets_to_remove.append(met)
            case "reaction":
                rxn = model.reactions.get_by_id(row["id"])
                rxns_to_remove.append(rxn)

    return mets_to_remove, rxns_to_remove

def update_gprs(model, gene_rules):
    for row in gene_rules.iter_rows(named=True):
        try:
            rxn = model.reactions.get_by_id(row["rxn_id"])
            rxn.gene_reaction_rule = row["new_gpr"]
        except KeyError:
            print(f"{row["rxn_id"]} not in model")

def clean_genes(model):
    genes_to_remove = [g.id for g in model.genes if not bool(g.reactions)]
    for gid in genes_to_remove:
        model.genes.remove(gid)

def update_model(model, updates_path):
    new_model = model.copy()
    metabolites = pl.read_excel(updates_path, sheet_name="metabolites")
    reactions = pl.read_excel(updates_path, sheet_name="reactions")
    boundary = pl.read_excel(updates_path, sheet_name="boundary_rxns")
    gene_rules = pl.read_excel(updates_path, sheet_name="gpr")
    deletions = pl.read_excel(updates_path, sheet_name="deletions")

    mets_to_add = update_metabolites(metabolites)
    rxns_to_add = update_reactions(reactions)
    mets_to_remove, rxns_to_remove = update_dels(new_model, deletions)

    new_model.remove_metabolites(mets_to_remove)
    new_model.remove_reactions(rxns_to_remove)

    print("Adding mets")
    new_model.add_metabolites(mets_to_add)
    print("Adding rxns")
    new_model.add_reactions(rxns_to_add)
    for row in reactions.iter_rows(named=True):
        rxn = new_model.reactions.get_by_id(row["rxn_id"])
        rxn.build_reaction_from_string(row["reaction"])
    print("Adding boundaries")
    update_boundary_rxns(new_model, boundary)
    print("Adding gprs")
    update_gprs(new_model, gene_rules)

    return new_model

def remove_orphan_metabolites(model):
    rxns_to_remove = True
    removed_rxns = []
    removed_mets = []
    while rxns_to_remove:
        # orphans, deadends = find_blocked_mets(model)

        rxns_to_remove = []
        for met in model.metabolites:
            # met = model.metabolites.get_by_id(mid)
            if len(met.reactions) <= 1:
                rxns_to_remove.extend([*met.reactions])

        rxns_to_remove = set(rxns_to_remove)
        model.remove_reactions(rxns_to_remove)
        removed_rxns.extend(rxns_to_remove)

        mets_to_remove = [met for met in model.metabolites if len(met.reactions)==0]
        model.remove_metabolites(mets_to_remove)
        removed_mets.extend(mets_to_remove)

def qc_filter_dl(dl_df):
    return (
        dl_df
        .unpivot(
            index=["Protein_ID", "Localizations", "Signals"],
            variable_name="location",
            value_name="confidence",
        )
        #.filter(pl.col("confidence")>=0.4)
        .with_columns(
            (pl.col("confidence").max() - pl.col("confidence")).over("Protein_ID").alias("max_dif"),
            pl.col("location").str.to_lowercase(),
        )
        .filter(pl.col("max_dif")<=0.08)
        .sort(["Protein_ID", "confidence"], descending=True)
        .with_columns(
            pl.col("location").replace_strict(compartments, default=None)
        )
        #.drop_nulls()
        .group_by("Protein_ID")
        .agg(pl.col("location"))
        .with_columns(
            pl.col("location").list.join(" or ")
        )
    )

def clean_gprs(model):
    reactions_list = []
    for rxn in model.reactions:
        new_dict = {
            "id": rxn.id,
            "gpr": rxn.gpr.to_string(),
        }
        reactions_list.append(new_dict)

    reactions = pl.DataFrame(reactions_list, orient="row")
    query = (
        reactions
        .with_columns(
            pl.col("gpr").str.split(" or ").alias("genes")
        )
        .with_columns(
            pl.col("genes").list.unique().alias("u_genes")
        )
        .with_columns(
            (pl.col("genes").list.len() != pl.col("u_genes").list.len()).alias("repeated")
        )
        .filter(pl.col("repeated"))
        .with_columns(
            pl.col("u_genes").list.join(" or ").alias("new_gpr")
        )
    )

    for row in query.iter_rows(named=True):
        rxn = model.reactions.get_by_id(row["id"])
        rxn.gene_reaction_rule = row["new_gpr"]

def get_reaction_annotation(reactions):
    query = (
        reactions
        .select(["id", "metanetx.reaction", "seed.reaction", "metacyc.reaction", "rhea", "reactome"])
    )
    rxns_wo_mtntx = query.filter(pl.col("metanetx.reaction").is_null())
    bad_mtntx = (
        query
        .filter(~pl.col("metanetx.reaction").is_null())
        .join(mtntx_db, left_on="metanetx.reaction", right_on="id", how="anti")
        .with_columns(pl.lit(True).alias("bad_mtntx"))
    )
    rxns_w_mtntx = (
        query.select(["id", "metanetx.reaction"])
        .filter(~pl.col("metanetx.reaction").is_null())
        .join(mtntx_db, left_on="metanetx.reaction", right_on="id", how="inner")
    )
    return pl.concat([rxns_wo_mtntx, bad_mtntx, rxns_w_mtntx], how="diagonal")

if __name__=="__main__":
    # Inputs
    model_path = Path("models/draft/v0.0.3/nies144/nies144.xml")
    excel_path = Path("models/draft/v0.0.3/nies144/nies144.xlsx")
    updates_path = Path("data/1_interim/curation/updates_to_model.xlsx")

    # Load Model
    base = cobra.io.read_sbml_model(model_path)
    reactions = pl.read_excel(excel_path, sheet_name="reactions")

    # Update the model
    model = update_model(base, updates_path)

    sol = model.optimize()
    model.summary(solution=sol)

    # Find and clean gprs with repeated genes
    clean_gprs(model)
    clean_genes(model)

    # Fix annotations
    # Fix reaction annotation
    # TO DO: move this to another step (maybe 2?)
    uris_to_fix = {
        "ECNumber": "ec-code",
        "KEGGReaction": "kegg.reaction",
        "biocyc": "metacyc.reaction",
        "ReactomeReaction": "reactome",
    }
    for rxn in model.reactions:
        for key, val in uris_to_fix.items():
            identifier = rxn.annotation.get(key, None)
            if identifier:
                rxn.annotation[val] = identifier
                rxn.annotation.pop(key)
    reactions, metabolites, genes = write_excel(model, "models/draft/v0.0.4/nies144/nies144.xlsx")
    annotation_df = get_reaction_annotation(reactions)

    for row in annotation_df.iter_rows(named=True):
        rxn_id = row["id"]
        row.pop("id")
        annotation = {key: val for key, val in row.items() if val is not None}
        rxn = model.reactions.get_by_id(rxn_id)
        rxn.annotation = annotation

    reactions, metabolites, genes = write_excel(model, "models/draft/v0.0.4/nies144/nies144.xlsx")

    # Exogenous genes
    exog = [g.id for g in model.genes if "K" in g.id]
    exog.sort()
    len(exog)

    # Test remove reaction
    rxn_id = "GLYDHD"
    rxn = model.reactions.get_by_id(rxn_id)
    model.remove_reactions([rxn])
    model.slim_optimize()

    model =
    model_path = Path("models/draft/v0.0.2/nies.xml")
    model_path = Path("models/draft/v0.0.1/hlacustris.xml")

    model = cobra.io.read_sbml_model(model_path)
    model.reactions.get_by_id("AGMIS")

    import os
    rxns_to_test = []
    with Path("test_rmv").open("r") as f:
        for line in f:
            rxns_to_test.append(line.strip())

    rxns_to_rmv = [base.reactions.get_by_id(rid) for rid in rxns_to_test]
    for rxn in rxns_to_rmv:
        base.remove_reactions([rxn])
        print(rxn.id, base.slim_optimize())
