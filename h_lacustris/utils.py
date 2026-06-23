from cobra import Metabolite, Reaction

import polars as pl
from xlsxwriter import Workbook


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

file = "data/external/bigg_metabolites.tsv"
bigg_metabolites = pl.read_csv(file, separator="\t", encoding="utf8-lossy", ignore_errors=True)
file = "data/external/bigg_reactions.tsv"
bigg_reactions = pl.read_csv(file, separator="\t", encoding="utf8-lossy", ignore_errors=True)

def update_metabolites(metabolites_df: pl.DataFrame) -> list[Metabolite]:
    """Create list of metabolite instances to add to a model."""
    mets_to_add = []
    for met in metabolites_df.iter_rows(named=True):
        # Get met info from the database
        new_met = Metabolite(
            id=met["met_id"],
            formula=met["formula"],
            name=met["name"],
            charge=met["charge"],
            compartment=met["met_id"][-1],
        )

        # Only the annotations should be left in the row dict
        keys_to_remove = ["met_id", "formula", "name", "charge", "met_id"]
        for key in keys_to_remove:
            met.pop(key)
        new_met.annotations = met

        mets_to_add.append(new_met)

    return mets_to_add

def update_reactions(reactions_df: pl.DataFrame) -> list[Reaction]:
    """Create list of metabolite instances to add to a model."""
    rxns_to_add = []
    for rxn in reactions_df.iter_rows(named=True):
        # Get reaction info from the dataframe
        new_rxn = Reaction(
            id=rxn["rxn_id"],
            name=rxn["name"],
            subsystem=rxn["subsystem"],
            lower_bound=rxn["lower_bound"],
            upper_bound=rxn["upper_bound"],
        )
        if rxn["gpr"]:
            new_rxn.gene_reaction_rule = rxn["gpr"]

        # Get bigg database_links
        annotation = bigg_reactions.filter(
            pl.col("bigg_id") == rxn["rxn_id"]
        )["database_links"].to_list()
        # A patch for db_links with null values
        if annotation and annotation[0] is None:
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
    return rxns_to_add

def write_excel(model: cobra.Model, file_name: str) -> None:
    reactions = []
    for rxn in model.reactions:
        this_rxn = {
            "id": rxn.id,
            "name": rxn.name,
            "reaction": rxn.build_reaction_string(),
            "gpr": rxn.gpr.to_string(),
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

    genes = []
    for gen in model.genes:
        this_gene = {
            "id": gen.id,
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

def update_gene_ids(model, gene_map):
    print(f"\nUpdating genes in {model.name}")
    for rxn in tqdm(model.reactions):
        rule = rxn.gene_reaction_rule
        if not rule:
            continue
        gpr = rxn.gpr
        for gene in gpr.genes:
            if not gene in gene_map:
                continue
            old_id = gene
            new_id = gene_map[old_id]
            rule = re.sub(rf"\b{re.escape(old_id)}\b", new_id, rule)

        rxn.gene_reaction_rule = rule
        rxn.update_genes_from_gpr()

    print(f"Total genes before: {len(model.genes)}")
    genes_to_remove = [gene.id for gene in model.genes
                               if not bool(gene.reactions)]
    print(f"Removing {len(genes_to_remove)} unused genes")
    for gene_id in genes_to_remove:
        model.genes.remove(gene_id)
    print(f"Total genes after update: {len(model.genes)}")
