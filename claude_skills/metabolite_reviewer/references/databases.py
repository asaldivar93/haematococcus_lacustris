from pathlib import Path
import polars as pl

OBSOLET_RXNS = ["MNXR01", "MNXR02", "MNXR03"]
DATABASE_PATH = Path("models/external/databases/")

# Load MetaNetX database
# Reaction Info
metanetx_rxns_path = DATABASE_PATH / "metanetx/reac_prop.tsv"
metanetx_rxns = pl.read_csv(
    metanetx_rxns_path,
    separator="\t",
    comment_prefix="#",
    has_header=False,
    new_columns=["id", "reaction", "xref", "ec-code", "balanced", "transport"],
)

# Reaction to external databases
metanetx_xref_path = DATABASE_PATH / "metanetx/reac_xref.tsv"
metanetx_xref = pl.read_csv(
    metanetx_xref_path,
    separator="\t",
    comment_prefix="#",
    has_header=False,
    ignore_errors=True,
    new_columns=["xref", "id", "equation"],
)

# Metabolite Info
metanetx_met_path = DATABASE_PATH / "metanetx/chem_prop.tsv"
metanetx_mets = pl.read_csv(
    metanetx_met_path,
    separator="\t",
    comment_prefix="#",
    has_header=False,
    ignore_errors=True,
    new_columns=[
        "id",
        "name",
        "xref",
        "formula",
        "charge",
        "mass",
        "inchi",
        "inchikey",
        "smiles",
    ],
)

# Metabolites to external databases
metanetx_met_xref_path = DATABASE_PATH / "metanetx/chem_xref.tsv"
metanetx_mets_xref = pl.read_csv(
    metanetx_met_xref_path,
    separator="\t",
    comment_prefix="#",
    quote_char=None,
    has_header=False,
    ignore_errors=True,
    new_columns=["xref", "id", "name"],
)

# Compartments used by MetaNetX
compartments_path = DATABASE_PATH / "metanetx/comp_prop.tsv"
comps_df = pl.read_csv(
    compartments_path,
    separator="\t",
    has_header=False,
    comment_prefix="#",
    ignore_errors=True,
    new_columns=["id", "name", "go"],
)
MTNTX_COMPARTMENTS = comps_df["name"].to_list()

# Transform metanetx_xref from long to wide format
# Translate xref names using identifiers.org notation
IDENTIFIERS_ORG_MAP = {
    "biggR": "bigg.reaction",
    "rhea": "rh",
    "seedR": "seed.reaction",
    "metacycR": "metacyc.reaction",
    "vmhR": "vmhreaction",
    "sabiorkR": "sabiork.reaction",
    "keggR": "kegg.reaction",
    "ReactomeReaction": "reactome",
}

obsolete_rxn = "secondary/obsolete/fantasy identifier"

rxn_xref = (
    metanetx_xref
    # Remove dummy reactions
    .filter(
        ~pl.col("equation").str.contains(obsolete_rxn),
        ~(pl.col("id") == "EMPTY"),
        ~pl.col("id").is_in(OBSOLET_RXNS),
    )
    # Replace xref ids with identifier.org uris
    .with_columns(
        pl.col("xref").str.replace_many(IDENTIFIERS_ORG_MAP),
    )
    # Get all external identifiers for each reaction
    .group_by("id")
    .agg(
        [
            pl.col("xref").unique().sort(),
        ],
    )
    .explode("xref")
    # Create a row with xref to MetaNetX id relationships
    .with_columns(
        pl.col("xref")
        .str.split_exact(":", 1)
        .struct.rename_fields(["database", "xid"]),
    )
    .unnest("xref")
    .group_by("id", "database")
    .agg(pl.col("xid"))
    .with_columns(pl.col("xid").list.join("|"))
    # Make the dataframe wide format
    .pivot(
        on="database",
        index="id",
        values="xid",
    )
)

METANETX_RXNS_DB = (
    metanetx_rxns
    .select("id", "ec-code")
    .join(rxn_xref, on="id", how="inner")
    .sort("id")
    .with_columns(pl.col("ec-code").str.replace(";", "|"))
    .rename({"rh": "rhea"})
)

# Transform metanetx_xref from long to wide format
# Translate xref names using identifiers.org notation
replace_mets = {
    "biggm": "bigg.metabolite",
    "metacycm": "metacyc.compound",
    "keggd": "kegg.drug",
    "keggc": "kegg.compound",
    "keggg": "kegg.glycan",
    "rheap": "rhea",
    "rheag": "rhea",
    "lipidmapsm": "lipidmaps",
    "seedm": "seed.compound",
    "vmhmetabolite": "vmhm",
    "sabiorkm": "sabiork.compound",
    "reactomem": "reactome",
}
rmv_mets = ["WATER", "MNXM01", "MNXM02", "MNXM03", "BIOMASS", "MNXM1"]
rmv_databases = ["envipathm", "envipath"]

mets_xref = (
    metanetx_mets_xref.filter(
        ~pl.col("id").is_in(rmv_mets),
    )
    # Flatten database names
    .with_columns(pl.col("xref").str.to_lowercase())
    # Replace xref with identifiers.org
    .with_columns(
        pl.col("xref").str.replace_many(replace_mets),
    )
    .group_by("id")
    .agg(
        [
            pl.col("xref").unique().sort(),
        ],
    )
    .explode("xref")
    .with_columns(
        pl.col("xref")
        .str.split_exact(":", 1)
        .struct.rename_fields(["database", "xid"]),
    )
    .unnest("xref")
    .filter(~pl.col("database").is_in(rmv_databases))
    .group_by("id", "database")
    .agg(pl.col("xid"))
    .with_columns(pl.col("xid").list.sort().list.join("|"))
    .pivot(
        on="database",
        index="id",
        values="xid",
    )
)

METANETX_METS_DB = (
    metanetx_mets
    .with_columns(pl.col("inchi").str.replace("InChI=", ""))
    .join(mets_xref, on="id", how="inner")
    .select(pl.exclude("xref", "mass"))
    .sort("id")
)

bigg_rxns = pl.read_csv(
    DATABASE_PATH / "bigg/bigg_models_reactions.txt",
    separator="\t",
    ignore_errors=True,
    new_columns=["id"],
).rename({"reaction_string": "reaction"})

BIGG_RXNS_DB = (
    bigg_rxns
    .with_columns(pl.col("database_links").str.split(";"))
    .explode("database_links")
    .with_columns(
        pl.col("database_links")
        .str.extract(r"(https?://.*)")
        .str.replace("META:","")
        .str.replace("biocyc", "metacyc.reaction")
        .str.split("/")
        .list.tail(2)
        .list.to_struct(fields=["database", "xref"]),
    )
    .unnest("database_links")
    .group_by("id", "database")
    .agg(pl.col("xref"))
    .with_columns(pl.col("xref").list.sort().list.join("|"))
    .pivot(
        index="id",
        on="database",
        values="xref",
    )
    .drop("null")
    .join(
        bigg_rxns.select(["id", "name", "reaction", "model_list", "old_bigg_ids"]),
        on="id",
        how="inner",
    )
    .with_columns(pl.col("id").alias("bigg.reaction"))
)

bigg_mets = pl.read_csv(
    DATABASE_PATH / "bigg/bigg_models_metabolites.txt",
    separator="\t",
    ignore_errors=True,
).rename({"universal_bigg_id": "id"})

BIGG_METS_DB = (
    bigg_mets
    .with_columns(pl.col("database_links").str.split(";"),)
    .unique("id")
    .explode("database_links")
    .with_columns(
        pl.col("database_links")
        .str.replace("biocyc", "metacyc.compound")
        .str.replace("META:","")
        .str.replace("CHEBI:","")
        .str.extract(r"(https?://.*)")
        .str.split("/")
        .list.tail(2)
        .list.to_struct(fields=["database", "xref"]),
    )
    .unnest("database_links")
    .group_by("id", "database")
    .agg(pl.col("xref"))
    .with_columns(pl.col("xref").list.sort().list.join("|"))
    .sort("id")
    .pivot(
        index="id",
        on="database",
        values="xref",
    )
    .drop("null")
    .join(
        bigg_mets.select("id", "name", "model_list", "old_bigg_ids"),
        on="id",
        how="inner"
    )
    .unique("id")
    .with_columns(pl.col("id").alias("bigg.metabolite"))
)

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
