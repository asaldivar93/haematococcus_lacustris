
import cobra
import polars as pl
from cobra import Metabolite, Model, Reaction
from xlsxwriter import Workbook

compartments = {
    "u": "thylakoid",
    "m": "mitochondrion",
    "h": "chloroplast",
    "x": "glyoxysome",
    "c": "cytoplasm",
    "e": "extracellular"
}

replace_map = {
    "[u]": "_u",
    "[m]": "_m",
    "[h]": "_h",
    "[x]": "_x",
    "[c]": "_c",
    "[e]": "_e",
}

replace_map_rxns = {
    "(e)": "_e",
    "(c)": "_c",
    "(h)": "_h",
    "(m)": "_m",
    "(x)": "_x",
    "(u)": "_u",
}

db_key_map = {
    "KEGGCompound": "kegg.compound",
    "CHEBI": "chebi",
    "InChIKey": "inchikey",
    "HumanMetabolomeDatabase": "hmdb",
    "LipidMaps": "lipidmaps",
    "BioCyc": "biocyc",
    "ReactomeCompound": "reactome",
    "MetaNetX(MNX)Chemical": "metanetx.chemical",
    "SEEDCompound": "seed.compound",
    "RHEA": "rhea",
    "MetaNetX(MNX)Equation": "metanetx.reaction",
    "SEEDReaction": "seed.reaction",
    "Reactome Reaction": "reactome",
    "EC Number": "ec-code",
}

illegal_chars = ["-", "(", ")"]

file = "data/external/bigg_metabolites.tsv"
bigg_metabolites = pl.read_csv(file, separator="\t")
unique_bigg_mets_df = bigg_metabolites.unique("universal_bigg_id", keep="first")

def correct_met_ids(mets_in_model_df):
    """Change met ids with match to bigg ids."""
    # Get the list of unique universal ids in bigg and in the model
    mets_in_model_df = metabolites_df
    unique_mets_df = mets_in_model_df.unique("universal_id", keep="first")
    # Find the metabolite ids not found in bigg
    mets_with_bigg = unique_mets_df["universal_id"].is_in(unique_bigg_mets_df["universal_bigg_id"])
    mets_wo_bigg = unique_mets_df.filter(~mets_with_bigg).select(
        pl.col("universal_id"),
        pl.col("name"),
    )

    # Find bigg metabolites by name
    mets_wo_bigg = mets_wo_bigg.with_columns(name=pl.col("name").str.to_lowercase())
    mets_with_name = mets_wo_bigg.join(
        unique_bigg_mets_df.with_columns(name=pl.col("name").str.to_lowercase()),
        on="name",
        how="inner",
    ).sort("bigg_id")
    has_zero_duplicates = mets_with_name.filter(
        mets_with_name["universal_id"].is_duplicated(),
    ).is_empty()
    # Verify there are not duplicates
    if not has_zero_duplicates:
        with pl.Config(tbl_rows=-1):
            print(
                mets_with_name.filter(
                    mets_with_name["universal_id"].is_duplicated(),
                ).sort("universal_id")
            )
        raise ValueError("Some mets have multiple matches to bigg ids")

    mets_wo_name = mets_wo_bigg.join(
        unique_bigg_mets_df.with_columns(name=pl.col("name").str.to_lowercase()),
        on="name",
        how="anti"
    )

    has_all_mets = mets_with_name.shape[0] + mets_wo_name.shape[0] == mets_wo_bigg.shape[0]
    if not has_all_mets:
        raise ValueError("Something is wrong")

    return mets_with_name, mets_wo_name

def correct_rxn_ids(rxns_in_model_df):
    rxns_in_bigg = rxns_in_model_df["id"].is_in(
        bigg_reactions["bigg_id"].append(new_rxn_ids["new_id"])
    )
    rxns_wo_bigg = rxns_in_model_df.filter(~rxns_in_bigg).select(
        pl.col("id"),
        pl.col("name"),
        pl.col("reaction"),
    )
    rxns_wo_bigg = rxns_wo_bigg.with_columns(name=pl.col("name").str.to_lowercase())
    rxns_with_name = rxns_wo_bigg.join(
        bigg_reactions.with_columns(name=pl.col("name").str.to_lowercase()),
        on="name",
        how="inner"
    )
    has_zero_duplicates = rxns_with_name.filter(
        rxns_with_name["id"].is_duplicated(),
    ).is_empty()
    if not has_zero_duplicates:
        with pl.Config(tbl_rows=-1):
            print(
                rxns_with_name.filter(
                    rxns_with_name["id"].is_duplicated(),
                ).sort("id")
            )
        raise ValueError("Some rxns have multiple matches to bigg ids")
    rxns_wo_name = rxns_wo_bigg.join(
        bigg_reactions.with_columns(name=pl.col("name").str.to_lowercase()),
        on="name",
        how="anti"
    )
    has_all_rxns = rxns_with_name.shape[0] + rxns_wo_name.shape[0] == rxns_wo_bigg.shape[0]
    if not has_all_rxns:
        raise ValueError("Something is wrong")

    return rxns_with_name, rxns_wo_name

def write_excel(model, file_name):
    reactions = []
    for rxn in model.reactions:
        this_rxn = {
            "id": rxn.id,
            "name": rxn.name,
            "reaction": rxn.build_reaction_string(),
            "gpr": rxn.gpr,
            "subsystem": rxn.subsystem,
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

    with Workbook(file_name) as wb:
        pl.DataFrame(reactions).write_excel(
            workbook=wb,
            worksheet="reactions",
        )
        pl.DataFrame(metabolites).write_excel(
            workbook=wb,
            worksheet="metabolites",
        )

if __name__=="__main":
    # Load bigg metabolites
    file = "data/external/bigg_metabolites.tsv"
    bigg_metabolites = pl.read_csv(file, separator="\t")
    # Load bigg reactions
    file = "data/external/bigg_reactions.tsv"
    bigg_reactions = pl.read_csv(file, separator="\t")
    # Load dict of manully curated metabolite ids
    file = "data/interim/curation/alexis/updates_to_model.xlsx"
    corrected_met_ids = pl.read_excel(file, sheet_name="corrected_ids").filter(
        pl.col("type") == "metabolite",
    )
    corrected_rxn_ids = pl.read_excel(file, sheet_name="corrected_ids").filter(
        pl.col("type") == "reaction",
    )
    # load dict of new_bigg_ids
    file = "data/interim/curation/alexis/updates_to_model.xlsx"
    new_met_ids = pl.read_excel(file, sheet_name="new_bigg_ids").filter(
        pl.col("type") == "metabolite",
    )
    new_rxn_ids = pl.read_excel(file, sheet_name="new_bigg_ids").filter(
        pl.col("type") == "reaction",
    )
    # append new ids to the replacement map for metabolites

    replace_map.update(dict(zip(corrected_met_ids["old_id"], corrected_met_ids["new_id"])))
    replace_map.update(dict(zip(new_met_ids["old_id"], new_met_ids["new_id"])))
    replace_map_rxns.update(dict(zip(corrected_rxn_ids["old_id"], corrected_rxn_ids["new_id"])))
    replace_map_rxns.update(dict(zip(new_rxn_ids["old_id"], new_rxn_ids["new_id"])))
    # Load metabolites in model and update metabolite ids
    file = "models/draft/v0.0.0/Hlacustris.xlsx"
    metabolites_df = pl.read_excel(file, sheet_name="metabolites")
    metabolites_df.with_columns(
        old_id=pl.col("id"),
        universal_id=pl.col("id").str.head(-3)
    )
    metabolites_df = metabolites_df.with_columns(
        id = pl.col("id").str.replace_many(replace_map),
    )
    metabolites_df = metabolites_df.with_columns(
        universal_id = pl.Series([met[:-2] for met in metabolites_df["id"]]),
    )
    # Find mets without bigg id (serch by id and by name)
    mets_with_name, mets_wo_name = correct_met_ids(metabolites_df)
    # Save results to data/interim
    file = "data/interim/curation/alexis/corrected_mets_automatic.csv"
    mets_with_name.select(
        pl.col("universal_id").alias("old_id"),
        pl.col("universal_bigg_id").alias("new_id")
    ).write_csv(file)
    file = "data/interim/curation/alexis/mets_wo_bigg.csv"
    mets_wo_name.write_csv(file)
    # Update metabolites_df
    replace_map.update(
        dict(zip(mets_with_name["universal_id"], mets_with_name["universal_bigg_id"]))
    )
    metabolites_df = metabolites_df.with_columns(
        id = pl.col("id").str.replace_many(replace_map),
        universal_id = pl.col("universal_id").str.replace_many(replace_map)
    )
    # Find metabolites with illegal characters in the id
    file = "data/interim/curation/alexis/illegal_mets.csv"
    urgent_mets = [met_id for met_id in metabolites_df["universal_id"]
                          if any(char in met_id for char in illegal_chars)]
    metabolites_df.filter(pl.col("universal_id").is_in(urgent_mets)).select(
        pl.col("universal_id"),
        pl.col("name"),
        pl.col("formula"),
    ).write_csv(file)

    # Load reactions in model and update metabolite ids
    file = "models/draft/v0.0.0/Hlacustris.xlsx"
    replace_map_rxns.update(dict(zip(corrected_rxn_ids["old_id"], corrected_rxn_ids["new_id"])))
    replace_map_rxns.update(dict(zip(new_rxn_ids["old_id"], new_rxn_ids["new_id"])))

    reactions_df = pl.read_excel(
        file,
        sheet_name="reactions",
        schema_overrides={"lower_bound": pl.Float64, "upper_bound": pl.Float64}
    )
    reactions_df = reactions_df.with_columns(
        reaction = pl.col("reaction").str.replace_many(replace_map),
        id = pl.col("id").str.replace_many(replace_map_rxns),
    )
    reactions_df = reactions_df.with_columns(
        id = pl.col("id").str.replace("ARGDCI", "ARGDI")
    )

    # Find reactions without bigg id (search by id and by name)
    rxns_with_name, rxns_wo_name = correct_rxn_ids(reactions_df)
    # Save results to data/interim
    file = "data/interim/curation/alexis/corrected_rxns_automatic.csv"
    rxns_with_name.select(
        pl.col("id").alias("old_id"),
        pl.col("bigg_id").alias("new_id")
    ).write_csv(file)
    file = "data/interim/curation/alexis/rxns_wo_bigg.csv"
    rxns_wo_name.write_csv(file)
    replace_map_rxns.update(
        dict(zip(rxns_with_name["id"], rxns_with_name["bigg_id"]))
    )
    reactions_df = reactions_df.with_columns(
        id=pl.col("id").str.replace_many(replace_map_rxns)
    )
    # Find reactions with illegal characters in the id
    file = "data/interim/curation/alexis/illegal_rxns.csv"
    urgent_rxns = [rxn_id for rxn_id in reactions_df["id"]
                          if any(char in rxn_id for char in illegal_chars)]
    reactions_df.filter(pl.col("id").is_in(urgent_rxns)).select(
        pl.col("id"),
        pl.col("name"),
        pl.col("reaction"),
    ).write_csv(file)

    # Filter gprs to review
    # Reactions with dunaliella genes without haem genes
    duns_gprs = reactions_df.filter(
        pl.col("Gene").is_null()
    ).filter(
        ~pl.col("Dunaliella Gene").is_null()
    )["id"].to_list()

    # Protein complexes
    complex_gprs = reactions_df.filter(
        pl.col("Dunaliella Gene").str.contains("and")
    )["id"].to_list()

    file = "data/interim/curation/alexis/gprs_to_check.csv"
    gprs_to_check = set(duns_gprs + complex_gprs)
    reactions_df.filter(pl.col("id").is_in(gprs_to_check)).join(
        bigg_reactions.rename({"bigg_id": "id"}),
        on="id",
        how="left"
    ).select(
        pl.col("id"),
        pl.col("name"),
        pl.col("reaction"),
        pl.col("Gene"),
        pl.col("Dunaliella Gene"),
        pl.col("ec-code"),
        pl.col("database_links")
    ).write_csv(file)

    # Create model and assing compartments
    model = Model("iHlct", "H. lacustris v0.0.1")
    model.compartments = compartments

    # Add metabolites to model
    mets_to_add = []
    for met in metabolites_df.iter_rows(named=True):
        new_met = Metabolite(
            id=met["id"],
            formula=met["formula"],
            name=met["name"],
            charge=met["charge"],
            compartment=met["id"][-1],
        )

        # Get bigg database_links
        annotation = bigg_metabolites.filter(
            pl.col("universal_bigg_id") == met["universal_id"]
        )["database_links"].to_list()
        # A patch for db_links with null values
        if annotation:
            if annotation[0] is None:
                annotation = []
        # Add annotation if available
        if annotation:
            # get nested list of [[db_name, db_entry]]
            annotation_list = [entry.split(":", 1) for entry
                                                   in annotation[0].replace(" ","").split(";")]
            # Remove the uri from the db_entry
            annotation_df = pl.DataFrame(
                annotation_list,
                schema=["db_name", "identifier"],
                orient="row"
            ).unique("db_name", keep="first").with_columns(
                db_name=pl.col("db_name").str.replace_many(db_key_map),
                identifier=pl.col("identifier").str.split("/").list.get(-1),
            )
            # Build a dict of {db_name: db_entry}
            annotation_dict = dict(zip(annotation_df["db_name"], annotation_df["identifier"]))
            new_met.annotation = annotation_dict
        mets_to_add.append(new_met)

    model.add_metabolites(mets_to_add)

    # Add reactions to model
    rxns_to_add = []
    for rxn in reactions_df.iter_rows(named=True):
        new_rxn = Reaction(
            id=rxn["id"],
            name=rxn["name"],
            subsystem=rxn["subsystem"],
            lower_bound=rxn["lower_bound"],
            upper_bound=rxn["upper_bound"],
        )
        if rxn["Gene"]:
            new_rxn.gene_reaction_rule = rxn["Gene"]

        # Get bigg database_links
        annotation = bigg_reactions.filter(
            pl.col("bigg_id") == rxn["id"]
        )["database_links"].to_list()
        # A patch for db_links with null values
        if annotation:
            if annotation[0] is None:
                annotation = []
        # Add annotation if available
        if annotation:
            # get nested list of [[db_name, db_entry]]
            annotation_list = [entry.split(":", 1) for entry
                                                   in annotation[0].replace(" ","").split(";")]
            # Remove the uri from the db_entry
            annotation_df = pl.DataFrame(
                annotation_list,
                schema=["db_name", "identifier"],
                orient="row"
            ).unique("db_name", keep="first").with_columns(
                db_name=pl.col("db_name").str.replace_many(db_key_map),
                identifier=pl.col("identifier").str.split("/").list.get(-1),
            )
            # Build a dict of {db_name: db_entry}
            annotation_dict = dict(zip(annotation_df["db_name"], annotation_df["identifier"]))
            new_rxn.annotation = annotation_dict
        rxns_to_add.append(new_rxn)
    model.add_reactions(rxns_to_add)

    # Add reactions to model
    for rxn in reactions_df.iter_rows(named=True):
        new_rxn = model.reactions.get_by_id(rxn["id"])
        new_rxn.build_reaction_from_string(rxn["reaction"])

    write_excel(model, "models/draft/v0.0.1/hlacustris.xlsx")
    cobra.io.write_sbml_model(model, "models/draft/v0.0.1/hlacustris.xlm")

        file = "models/draft/v0.0.0/Hlacustris.xlsx"
        reactions_df = pl.read_excel(
            file,
            sheet_name="reactions",
            schema_overrides={"lower_bound": pl.Float64, "upper_bound": pl.Float64}
        )
        reactions_df = reactions_df.with_columns(
            reaction = pl.col("reaction").str.replace_many(replace_map),
            id = pl.col("id").str.replace_many(replace_map_rxns),
        )

# Load metabolites in model and update metabolite ids
file = "models/draft/v0.0.0/Hlacustris.xlsx"
m = pl.read_excel(file, sheet_name="metabolites")

m.filter(
    pl.col("name")=="alpha-D-Galactose 1-phosphate"
)

m.filter(
    pl.col("name")=="gal1p-L[c]"
)


metabolites_df.filter(
    pl.col("id")=="gal1p_c"
)

mask = reactions_df["id"].is_duplicated()
reactions_df.filter(mask).sort("id")[0, 1]

metabolites_df.filter(
    pl.col("id").str.contains("glyc")
)
