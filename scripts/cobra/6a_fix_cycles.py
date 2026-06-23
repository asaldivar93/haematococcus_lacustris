from pathlib import Path

import cobra
import polars as pl
# Load Reference
ref_path = Path("models/draft/v0.0.2/hlacustris.xml")
reference = cobra.io.read_sbml_model(ref_path)
ref_rxns = {rxn.id for rxn in reference.reactions}

# Load the model
model_path = Path("models/draft/v0.0.6/nies144/nies144.xml")
model = cobra.io.read_sbml_model(model_path)
# model_rxns = {rxn.id for rxn in model.reactions}
# rxns_not_in_ref = model_rxns - ref_rxns
# rxns_to_remove = [model.reactions.get_by_id(id) for id in rxns_not_in_ref]
# model.remove_reactions(rxns_to_remove)

rxns_to_remove = ["2OXOADOXm"]
model.remove_reactions(rxns_to_remove)

mets_to_remove = [met for met in model.metabolites if len(met.reactions) < 1]
model.remove_metabolites(mets_to_remove)

model.slim_optimize()

ref = "models/external/iLB1027_lipid.xml"
model = cobra.io.read_sbml_model(ref)
model.genes

table = pl.read_csv(
    "data/1_interim/genomes/dict_nies_to_red.csv",
    separator = "\t",
    has_header = False,
    new_columns = ["query", "target", "pident", "evalue", "bits","qcov","tcov","qlen","tlen", "alln"]
)
cdd = pl.read_csv(
    "data/1_interim/genomes/cdd_nies_hits.tsv",
    separator = "\t",
    has_header = False,
    new_columns = ["query", "target", "pident", "evalue", "bits","qcov","tcov","qlen","tlen", "alln"]
)

gene = "KAJ9516274.1"
query = (
    table
    .filter(pl.col("query")==gene)
    .sort(["tlen", "pident"], descending=True)
)


domain = "PLN02344"
query2 = (
    cdd
    .filter(
        pl.col("query").is_in(query["target"].to_list()),
        pl.col("target")==domain,
    )
)

query.filter(pl.col("target")=="GFH24849.1")
