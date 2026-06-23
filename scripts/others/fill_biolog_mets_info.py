
import cobra
import polars as pl

from h_lacustris.databases import BIGG_METS_DB, BIGG_RXNS_DB

BIGG_METS_DB.collect()["database"].unique()
BIGG_RXNS_DB.collect()

BIGG_METS_DB.collect()

biolog_map = pl.read_excel("data/external/biolog_map.xlsx")
biolog_to_bigg = pl.read_excel("data/external/biolog_to_bigg_dict.xlsx")

query = (
    biolog_map
    .join(biolog_to_bigg, on=["plate", "metabolite"], how="left")
)
query.write_csv("data/external/biolog_map.csv")
