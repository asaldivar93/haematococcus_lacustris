from pathlib import Path

import polars as pl

OBSOLET_RXNS = ["MNXR01", "MNXR02", "MNXR03"]
DATABASE_PATH = Path("models/external/databases/")

# Load MetaNetX database
# Reaction Info
metanetx_rxns_path = DATABASE_PATH / "metanetx/reac_prop.tsv"
metanetx_rxns = pl.scan_csv(
    metanetx_rxns_path,
    separator="\t",
    has_header=False,
    new_columns=["id", "reaction", "xref", "ec-code", "balanced", "transport"],
)

# Reaction to external databases
metanetx_xref_path = DATABASE_PATH / "metanetx/reac_xref.tsv"
metanetx_xref = pl.scan_csv(
    metanetx_xref_path,
    separator="\t",
    has_header=False,
    ignore_errors=True,
    new_columns=["xref", "id", "equation"],
)

# Metabolite Info
metanetx_met_path = DATABASE_PATH / "metanetx/chem_prop.tsv"
metanetx_mets = pl.scan_csv(
    metanetx_met_path,
    separator="\t",
    has_header=False,
    ignore_errors=True,
    new_columns=[
        "id",
        "name",
        "xref",
        "formula",
        "charge",
        "mass",
        "inchl",
        "inchlkey",
        "smiles",
    ],
)

# Metabolites to external databases
metanetx_met_xref_path = DATABASE_PATH / "metanetx/chem_xref.tsv"
metanetx_mets_xref = pl.scan_csv(
    metanetx_met_xref_path,
    separator="\t",
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
}

obsolete_rxn = "secondary/obsolete/fantasy identifier"
xref_lazy = (
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
    # Make the dataframe wide format
    .pivot(
        on="database",
        on_columns=IDENTIFIERS_ORG_MAP.values(),
        index="id",
        values="xid",
        aggregate_function="first",
    )
)

METANETX_RXNS_DB = (
    metanetx_rxns.join(xref_lazy, on="id", how="inner")
    .sort("id")
    .with_columns(pl.col("ec-code").str.split(";"))
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

xref_mets_lazy = (
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
    .sort("id")
    .filter(
        ~pl.col("database").is_in(rmv_databases),
    )
    .pivot(
        on="database",
        on_columns=list(set(replace_mets.values())),
        index="id",
        values="xid",
        aggregate_function="first",
    )
)

METANETX_METS_DB = metanetx_mets.join(
    xref_mets_lazy,
    on="id",
    how="inner",
).sort("id")

bigg_lazy = pl.scan_csv(
    DATABASE_PATH / "bigg/bigg_reactions.tsv",
    separator="\t",
    ignore_errors=True,
    new_columns=["bigg.reaction"],
)

BIGG_RXNS_DB = bigg_lazy.select("bigg.reaction", "name", "reaction_string")

bigg_mets = pl.scan_csv(
    DATABASE_PATH / "bigg/bigg_metabolites.tsv",
    separator="\t",
    ignore_errors=True,
)
on_columns = list(set(replace_mets.values()))
on_columns.append("metanetx.chemical")
BIGG_METS_DB = (
    bigg_mets.with_columns(
        pl.col("database_links").str.split(";"),
    )
    .explode("database_links")
    .with_columns(
        pl.col("database_links")
        .str.replace("biocyc", "metacyc.compound")
        .str.extract(r"(https?://.*)")
        .str.split("/")
        .list.tail(2)
        .list.to_struct(fields=["database", "xref"]),
    )
    .unnest("database_links")
    # .filter(~pl.col("database").is_null())
    .pivot(
        on="database",
        on_columns=["metanetx.chemical"],
        values="xref",
        aggregate_function="first",
    )
)
