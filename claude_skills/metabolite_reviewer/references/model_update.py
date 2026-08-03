
import cobra
import polars as pl

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

def clean_metabolites(model):
    """Remove metabolites without reactions."""
    mets_to_remove = [met for met in model.metabolites if len(met.reactions)==0]
    model.remove_metabolites(mets_to_remove)
    # Remove orphaned reactions
    rxns_to_remove = [rxn for rxn in model.reactions
        if len(rxn.metabolites)==1 and not rxn.boundary]
    model.remove_reactions(rxns_to_remove)

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

def clean_genes(model):
    genes_to_remove = [g.id for g in model.genes if not bool(g.reactions)]
    for gid in genes_to_remove:
        model.genes.remove(gid)
