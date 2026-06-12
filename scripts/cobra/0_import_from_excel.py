"""Second try :'(."""

import cobra
import polars as pl
from cobra import Metabolite, Model, Reaction
from labutils.cobra.utils import cobra_to_excel

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

def replace_ids(id_dataframe, replace_map):
    unique_id = list(id_dataframe["unique_id"])

    for i, met in enumerate(unique_id):
        if met in replace_map:
            unique_id[i] = replace_map[met]

    return id_dataframe.with_columns(
        unique_id=pl.Series(unique_id)
    )

def find_ids_wo_bigg(id_dataframe):
    # replace old ids with new new ids
    # I'm using this method instead of replace_many to make sure that
    # replacements are an exact match, replace_many changes substrings
    # and was replacing ids in unexpected ways
    unique_ids_df = id_dataframe.unique("unique_id", keep="first")
    ids_with_bigg = unique_ids_df["unique_id"].is_in(unique_bigg_mets_df["universal_bigg_id"])
    mets_wo_bigg = unique_ids_df.filter(~ids_with_bigg).select(
        pl.col("unique_id"),
        pl.col("name"),
    )
    mets_wo_bigg = mets_wo_bigg.with_columns(name=pl.col("name").str.to_lowercase())
    mets_with_name = mets_wo_bigg.join(
        unique_bigg_mets_df.with_columns(name=pl.col("name").str.to_lowercase()),
        on="name",
        how="inner",
    )
    has_duplicates = sum(mets_with_name["unique_id"].is_duplicated())
    if has_duplicates:
        with pl.Config(tbl_rows=-1):
            print(
                mets_with_name.filter(
                    mets_with_name["unique_id"].is_duplicated(),
                ).sort("unique_id")
            )
        raise ValueError("Some mets have multiple matches to bigg ids")


    return mets_wo_bigg, mets_with_name

def fix_metabolites(file):
    """Replace illegal character in metabolit ids.

    Replace with manually curated metabolite ids where available, or replace
    compartment with illegal chars to sbml compliant (ie. [c] to _c)
    """
    metabolites_df = pl.read_excel(file, sheet_name="metabolites")
    # Load dict of manully curated metabolite ids
    file = "data/interim/curation/alexis/updates_to_model.xlsx"
    # corrected_ids
    corrected_met_ids = pl.read_excel(file, sheet_name="corrected_ids").filter(
        pl.col("type") == "metabolite",
    )
    # new bigg ids
    new_bigg_ids = pl.read_excel(file, sheet_name="new_bigg_ids").filter(
        pl.col("type") == "metabolite",
    )
    # append both lists
    old_met_ids = list(corrected_met_ids["old_id"])
    b = list(new_bigg_ids["old_id"])
    old_met_ids.extend(b)

    new_met_ids = list(corrected_met_ids["new_id"])
    b = list(new_bigg_ids["new_id"])
    new_met_ids.extend(b)

    # Get the unique bigg id
    metabolites_df = metabolites_df.with_columns(
        compartment=pl.col("id").str.tail(3),
        unique_id=pl.col("id").str.head(-3),
    )

    # Replace manually corrected ids
    replace_mets = dict(zip(old_met_ids, new_met_ids))
    metabolites_df = replace_ids(metabolites_df, replace_mets)

    # Mets with bigg id identified from description
    mets_wo_bigg, mets_with_name = find_ids_wo_bigg(metabolites_df)
    mets_with_name.write_csv("data/interim/curation/corrected_mets_by_name.csv")

    # Replace ids based on bigg database
    # I actually don't remember why this is here
    replace_mets = dict(zip(mets_with_name["unique_id"], mets_with_name["universal_bigg_id"]))
    metabolites_df = replace_ids(metabolites_df, replace_mets)

    # Mets without bigg after looking through bigg database
    mets_wo_bigg, mets_with_name = find_ids_wo_bigg(metabolites_df)
    mets_wo_bigg.write_csv("data/interim/curation/alexis/mets_wo_bigg.csv")

    # These are metabolites with illegal characters
    urgent_mets = [met_id for met_id in metabolites_df["unique_id"]
                          if any(char in met_id for char in illegal_chars)]
    metabolites_df.filter(pl.col("unique_id").is_in(urgent_mets)).select(
        pl.col("unique_id"),
        pl.col("name"),
        pl.col("formula"),
    ).write_csv("data/interim/curation/alexis/illegal_mets.csv")

    return metabolites_df

def find_rxns_wo_bigg(id_dataframe):

    unique_ids_df = id_dataframe.unique("unique_id", keep="first")
    ids_with_bigg = unique_ids_df["unique_id"].is_in(bigg_reactions["bigg_id"])
    mets_wo_bigg = unique_ids_df.filter(~ids_with_bigg).select(
        pl.col("unique_id"),
        pl.col("name"),
    )
    mets_wo_bigg = mets_wo_bigg.with_columns(name=pl.col("name").str.to_lowercase())
    mets_with_name = mets_wo_bigg.join(
        bigg_reactions.with_columns(name=pl.col("name").str.to_lowercase()),
        on="name",
        how="inner",
    )

    return mets_wo_bigg, mets_with_name

def fix_reactions(file):
    """Replace illegal character in reaction ids.

    Replace with manually curated metabolite ids where available, or replace
    compartment with illegal chars to sbml compliant (ie. [c] to _c)
    """
    metabolites_df = pl.read_excel(file, sheet_name="reactions")
    # Load dict of manully curated metabolite ids
    file = "data/interim/curation/alexis/updates_to_model.xlsx"
    # corrected_ids
    corrected_met_ids = pl.read_excel(file, sheet_name="corrected_ids").filter(
        pl.col("type") == "reaction",
    )
    # new bigg ids
    new_bigg_ids = pl.read_excel(file, sheet_name="new_bigg_ids").filter(
        pl.col("type") == "reaction",
    )
    # append both lists
    old_met_ids = list(corrected_met_ids["old_id"])
    b = list(new_bigg_ids["old_id"])
    old_met_ids.extend(b)

    new_met_ids = list(corrected_met_ids["new_id"])
    b = list(new_bigg_ids["new_id"])
    new_met_ids.extend(b)

    # Get the unique bigg id
    metabolites_df = metabolites_df.with_columns(
        unique_id=pl.col("id")
    )

    replace_mets = dict(zip(old_met_ids, new_met_ids))
    metabolites_df = replace_ids(metabolites_df, replace_mets)

    mets_wo_bigg, mets_with_name = find_rxns_wo_bigg(metabolites_df)
    mets_with_name = mets_with_name.filter(~pl.col("unique_id").is_in(new_met_ids))
    has_duplicates = sum(mets_with_name["unique_id"].is_duplicated())
    if has_duplicates:
        with pl.Config(tbl_rows=-1):
            print(
                mets_with_name.filter(
                    mets_with_name["unique_id"].is_duplicated(),
                ).sort("unique_id")
            )
        raise ValueError("Some mets have multiple matches to bigg ids")
    mets_with_name.write_csv("data/interim/curation/alexis/corrected_rxns_by_name.csv")

    # replace old ids with new new ids
    # I'm using this method instead of replace_many to make sure that
    # replacements are an exact match, replace_many changes substrings
    # and was replacing ids in unexpected ways
    replace_mets = dict(zip(mets_with_name["unique_id"], mets_with_name["bigg_id"]))
    metabolites_df = replace_ids(metabolites_df, replace_mets)

    mets_wo_bigg, mets_with_name = find_rxns_wo_bigg(metabolites_df)
    mets_wo_bigg.write_csv("data/interim/curation/alexis/rxns_wo_bigg.csv")

    urgent_mets = [met_id for met_id in metabolites_df["unique_id"]
                          if any(char in met_id for char in illegal_chars)]
    metabolites_df.filter(pl.col("unique_id").is_in(urgent_mets)).select(
        pl.col("unique_id"),
        pl.col("name"),
    ).write_csv("data/interim/curation/alexis/illegal_rxns.csv")

    return metabolites_df

if __name__=="__main__":
    file = "data/external/bigg_metabolites.tsv"
    bigg_metabolites = pl.read_csv(file, separator="\t")
    unique_bigg_mets_df = bigg_metabolites.unique("universal_bigg_id", keep="first")
    # Load bigg reactions
    file = "data/external/bigg_reactions.tsv"
    bigg_reactions = pl.read_csv(file, separator="\t")

    file = "models/draft/v0.0.0/Hlacustris.xlsx"
    metabolites_df = fix_metabolites(file)
    reactions_df = fix_reactions(file)

    metabolites_df = metabolites_df.with_columns(
        compartment=pl.col("compartment").str.replace("[", "_", literal=True)
    )
    metabolites_df = metabolites_df.with_columns(
        compartment=pl.col("compartment").str.replace("]", "", literal=True)
    )
    metabolites_df = metabolites_df.with_columns(
        new_id=pl.concat_str(["unique_id", "compartment"])
    )
    replace_rxns = dict(zip(metabolites_df["id"], metabolites_df["new_id"]))

    reactions_df = reactions_df.with_columns(
        reaction=pl.col("reaction").str.replace_many(replace_rxns)
    )
    reactions_df = reactions_df.with_columns(
        lower_bound=pl.col("lower_bound").cast(pl.Float64),
        upper_bound=pl.col("upper_bound").cast(pl.Float64),
        Objective=pl.col("Objective").cast(pl.Float64),
    )
    metabolites_df = metabolites_df.unique(subset=["new_id"])

    # create the model
    model = Model("iHlct", "H. lacustris v0.0.1")
    model.compartments = compartments

    # Add metabolites to model
    mets_to_add = []
    for met in metabolites_df.iter_rows(named=True):
        new_met = Metabolite(
            id=met["new_id"],
            formula=met["formula"],
            name=met["name"],
            charge=met["charge"],
            compartment=met["new_id"][-1],
        )

        # Get bigg database_links
        annotation = bigg_metabolites.filter(
            pl.col("universal_bigg_id") == met["unique_id"]
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
    for rxn in reactions_df.iter_rows(named=True):
        if any([isinstance(rxn["lower_bound"], str), isinstance(rxn["upper_bound"], str)]):
            print(rxn["id"])
        new_rxn = Reaction(
            id=rxn["unique_id"],
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
        new_rxn = model.reactions.get_by_id(rxn["unique_id"])
        new_rxn.build_reaction_from_string(rxn["reaction"])

    model.objective = "BIOMASS_hlacus_auto"
    cobra.io.save_json_model(model, "models/draft/v0.0.1/hlacustris.json")

    model = cobra.io.load_json_model("models/draft/v0.0.1/hlacustris.json")
    cobra.io.write_sbml_model(model, "models/draft/v0.0.1/hlacustris.xml")
    cobra.io.validate_sbml_model("models/draft/v0.0.1/hlacustris.xml")
    model = cobra.io.read_sbml_model("models/draft/v0.0.1/hlacustris.xml")
    cobra_to_excel(model, "models/draft/v0.0.1/hlacustris.xlsx")
    model.optimize()
