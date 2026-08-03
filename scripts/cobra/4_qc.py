import asyncio

from pathlib import Path

import cobra
import polars as pl
from labutils.cobra.io import write_excel

from h_lacustris.databases import (
    BIGG_METS_DB,
    BIGG_RXNS_DB,
    METANETX_METS_DB,
    METANETX_RXNS_DB,
)
from h_lacustris.unichem import check_inchikeys_unichem_bulk

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
    annotation = pl.read_excel(updates_path, sheet_name="annotation")

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
    print("Adding new annotations")
    for row in annotation.iter_rows(named=True):
        try:
            match row["type"]:
                case "reaction":
                    item = new_model.reactions.get_by_id(row["id"])
                case "metabolite":
                    item = new_model.metabolites.get_by_id(row["id"])
                case _:
                    print(row)

            item.annotation[row["database"]] = row["xref"]
        except KeyError:
            print(f"{row["id"]} not in model")

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

def get_reaction_annotation(reactions_df):
    """Find rection in databases.

    Finds the reaction in bigg and adds the annotations from bigg. Then finds
    the reactions in metanetx and adds the annotations. Bigg annotations are
    superseded by metanetx.
    """
    new_annotations = [col for col in BIGG_RXNS_DB.columns
                       if col not in reactions_df.columns]
    query = (
        reactions_df
        .update(
            BIGG_RXNS_DB.select(pl.exclude(["reaction", "metanetx.reaction"])),
            on="id",
            how="left",
            include_nulls=False,
        )
        .join(
            BIGG_RXNS_DB.select("id", *new_annotations),
            on="id",
            how="left",
        )
        .drop("model_list", "old_bigg_ids")
    )
    new_annotations = [col for col in METANETX_RXNS_DB.columns
                       if col not in query.columns]
    return (
        query
        .update(
            METANETX_RXNS_DB
            .rename({"id": "metanetx.reaction"})
            .select(pl.exclude("bigg.reaction")),
            on="metanetx.reaction",
            how="left",
            include_nulls=False,
        )
        .join(
            METANETX_RXNS_DB.select("id", *new_annotations),
            left_on="metanetx.reaction",
            right_on="id",
            how="left",
        )
    )

def get_metabolite_annotation(metabolites_df):
    new_annotations = [col for col in BIGG_METS_DB.columns
                       if col not in metabolites_df.columns]
    query = (
        metabolites_df
        .rename({"id": "this_id"})
        .with_columns(pl.col("this_id").str.head(-2).alias("id"))
        .update(
            BIGG_METS_DB.select(pl.exclude("metanetx.chemical")),
            on="id",
            how="left",
            include_nulls=False,
        )
        .join(
            BIGG_METS_DB.select("id", *new_annotations),
            on="id",
            how="left",
        )
        .drop("id", "model_list", "old_bigg_ids")
        .rename({"this_id": "id"})
    )

    new_annotations = [col for col in METANETX_METS_DB.columns
                       if col not in query.columns]
    return (
        query
        .update(
            METANETX_METS_DB
            .rename({"id": "metanetx.chemical"})
            .select(pl.exclude("bigg.metabolite")),
            on="metanetx.chemical",
            how="left",
            include_nulls=False
        )
        .join(
            METANETX_METS_DB.select("id", *new_annotations),
            left_on="metanetx.chemical",
            right_on="id",
            how="left"
        )
    )

def update_annotations(model, metabolites_df, reactions_df):
    for row in metabolites_df.iter_rows(named=True):
        met = model.metabolites.get_by_id(row["id"])
        met.name = row["name"]
        met.formula = row["formula"]
        met.charge = row["charge"]
        for field in ["id", "name", "formula", "charge"]:
            row.pop(field)

        annotation = {key: value for key, value in row.items() if value is not None}
        met.annotation = annotation

    for row in reactions_df.iter_rows(named=True):
        rxn = model.reactions.get_by_id(row["id"])
        for field in ["id", "name", "reaction", "gpr", "subsystem"]:
            row.pop(field)
        annotation = {key: value for key, value in row.items() if value is not None}
        rxn.annotation = annotation

def clean_metabolites(model):
    """Remove metabolites without reactions."""
    mets_to_remove = [met for met in model.metabolites if len(met.reactions)==0]
    model.remove_metabolites(mets_to_remove)
    # Remove orphaned reactions
    rxns_to_remove = [rxn for rxn in model.reactions
        if len(rxn.metabolites)==1 and not rxn.boundary]
    model.remove_reactions(rxns_to_remove)

async def validate_inchikey(metabolites_df):
    inchi_list = metabolites_df["inchikey"].drop_nulls().unique().to_list()
    valid_inchi_dict = await check_inchikeys_unichem_bulk(inchi_list, checkpoint_path="inchi_checkpoint.json")
    valid_inchi_df = pl.DataFrame(valid_inchi_dict)
    return (
        metabolites_df
        .join(
            valid_inchi_df,
            on="inchikey",
            how="left"
        )
    )

def validate_metanetx(reactions_df, metabolites_df):
    query = (
        reactions_df
        .filter(~pl.col("metanetx.reaction").is_null())
        .select("id", "metanetx.reaction")
        .join(
            METANETX_RXNS_DB,
            left_on="metanetx.reaction",
            right_on="id",
            how="anti",
        )
        .with_columns(pl.lit(False).alias("valid_metanetx"))
        .drop("metanetx.reaction")
    )
    reactions_df = reactions_df.join(query, on="id", how="left")

    query = (
        metabolites_df
        .filter(~pl.col("metanetx.chemical").is_null())
        .select("id", "metanetx.chemical")
        .join(
            METANETX_METS_DB,
            left_on="metanetx.chemical",
            right_on="id",
            how="anti",
        )
        .with_columns(pl.lit(False).alias("valid_metanetx"))
        .drop("metanetx.chemical")
    )

    metabolites_df = metabolites_df.join(query, on="id", how="left")

    return reactions_df, metabolites_df

if __name__=="__main__":
    # Inputs
    model_path = Path("models/draft/v0.0.3/nies144/nies144.xml")
    excel_path = Path("models/draft/v0.0.3/nies144/nies144.xlsx")
    updates_path = Path("data/1_interim/curation/updates_to_model.xlsx")

    # Load Model
    base = cobra.io.read_sbml_model(model_path)

    # Update the model
    model = update_model(base, updates_path)

    sol = model.optimize()
    model.summary(solution=sol)

    # Clear dangling genes and metabolites
    clean_gprs(model)
    clean_genes(model)
    clean_metabolites(model)

    # Fix annotations
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

    uris_to_fix = {
        "biocyc": "metacyc.compound",
        "KEGGDrug": "kegg.drug",
        "KEGGGlycan": "keeg.glycan",
        "inchlkey": "inchikey"
    }
    for met in model.metabolites:
        for key, val in uris_to_fix.items():
            identifier = met.annotation.get(key, None)
            if identifier:
                met.annotation[val] = identifier
                met.annotation.pop(key)

    # update annotations
    out_path = "models/draft/v0.0.4/nies144/nies144.xlsx"
    reactions_df, metabolites_df = write_excel(model, out_path)
    reactions_df = get_reaction_annotation(reactions_df)
    metabolites_df = get_metabolite_annotation(metabolites_df)
    metabolites_df = await validate_inchikey(metabolites_df)
    reactions_df, metabolites_df = validate_metanetx(reactions_df, metabolites_df)

    # save the changes
    update_annotations(model, metabolites_df, reactions_df)
    reactions_df, metabolites_df  = write_excel(model, out_path)
    cobra.io.write_sbml_model(model, "models/draft/v0.0.4/nies144/nies144.xml")


    # Exogenous genes
    exog = [g.id for g in model.genes if "K" in g.id]
    exog.sort()
    len(exog)
    # Test remove reaction
    model = cobra.io.read_sbml_model("models/draft/v0.0.4/nies144/nies144.xml")
    test = model.copy()
    rxn_id = "NOS2"
    rxn = test.reactions.get_by_id(rxn_id)
    test.remove_reactions([rxn])
    sol = test.optimize()
    sol
    met = test.metabolites.get_by_id("man_c")
    test.add_boundary(met, type="sink")
    met.summary(solution=sol)

    model_path = Path("models/draft/v0.0.2/nies.xml")
    model_path = Path("models/draft/v0.0.1/hlacustris.xml")

    model = cobra.io.read_sbml_model(model_path)
    model.reactions.get_by_id("AGMIS")

    rxns_to_test = []
    with Path("test_rmv").open("r") as f:
        for line in f:
            rxns_to_test.append(line.strip())

    rxns_to_rmv = [test.reactions.get_by_id(rid) for rid in rxns_to_test]
    for rxn in rxns_to_rmv:
        test.remove_reactions([rxn])
        sol = test.slim_optimize()
        if sol <= 1e-6:
            test.add_reactions([rxn])
            print(rxn.id, sol)



    for row in annotation.iter_rows(named=True):
        match row["type"]:
            case "reaction":
                item = model.reactions.get_by_id(row["id"])
            case "metabolite":
                item = model.metabolites.get_by_id(row["id"])
            case _:
                print(row)
        #
        # for f in ["type", "id", "method"]:
        #     row.pop(f)
        #
        # for db, val in row.items():
        #     if val is not None:
        #         item.annotation[db] = val


row
