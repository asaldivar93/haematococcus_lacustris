from pathlib import Path

import cobra
import polars as pl
from labutils.cobra.io import write_excel

def clean_gprs(model):
    genes_to_remove = [g.id for g in model.genes if not bool(g.reactions)]
    for gid in genes_to_remove:
        model.genes.remove(gid)

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



if __name__=="__main__":
    model_path = Path("models/draft/v0.0.7/nies144/nies144.xml")
    excel_path = Path("models/draft/v0.0.7/nies144/nies144.xlsx")
    dl_nies_path = Path("data/2_processed/annotation/dl_nies144.csv")
    dl_red_path = Path("data/2_processed/annotation/dl_redball.csv")

    model = cobra.io.read_sbml_model(model_path)
    compartments = {item: key for key, item in model.compartments.items()}
    reactions = pl.read_excel(
        excel_path,
        sheet_name="reactions"
    )
    dl_nies = pl.read_csv(dl_nies_path)
    q = (
        dl_nies
        .with_columns()
    )

    dl_nies = qc_filter_dl(dl_nies)
    dl_nies_dict = dict(zip(dl_nies["Protein_ID"], dl_nies["location"]))
    dl_red = pl.read_csv(dl_red_path)
    dl_red = qc_filter_dl(dl_red).to_dict()
    dl_red_dict = dict(zip(dl_red["Protein_ID"], dl_red["location"]))

    # Find gprs with repeated genes
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
    clean_gprs(model)

    write_excel(model, "models/draft/v0.0.7/nies144/nies144.xlsx")

    reactions = pl.read_excel(
        excel_path,
        sheet_name="reactions"
    )

    # Add a list with compartments to the reaction
    compartments_by_rxn = []
    for rxn in model.reactions:
        new_dict = {
            "id": rxn.id,
            "compartments": list(rxn.compartments)
        }
        compartments_by_rxn.append(new_dict)
    reactions = (
        reactions
        .join(
            pl.DataFrame(compartments_by_rxn, orient="row"),
            on="id",
            how="inner"
        )
        .drop_nulls("gpr")
        .with_columns(
            pl.col("gpr").str.contains("KAJ").alias("exogenous")
        )
    )

    query = (
        reactions
        .with_columns(
            pl.col("gpr")
             .str.replace_many(dl_nies_dict)
             .str.replace_many(dl_red_dict)
             .alias("dl_location"),
        )
        .with_columns(
            pl.col("dl_location").str.split(" or ").list.unique()
        )
    )

    rows_to_append = []
    for row in query.iter_rows(named=True):
        has_loc = all([g in compartments.values() for g in row["dl_location"]])
        consitent = all([g in row["compartments"] for g in row["dl_location"]])
        row["all_loc"] = has_loc
        row["dl_consistent"] = consitent
        rows_to_append.append(row)

    query = pl.DataFrame(rows_to_append, orient="row").sort("id")
    query.write_excel("test.xlsx")
