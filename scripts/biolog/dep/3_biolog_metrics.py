import json
from pathlib import Path
from copy import deepcopy

import cobra
import polars as pl
#from scripts.utils import update_metabolites, update_reactions, write_excel
from memote.suite.tests.test_biomass import test_biomass_consistency
from sklearn.metrics import confusion_matrix, matthews_corrcoef
from xlsxwriter import Workbook


BIOLOG_TO_BIGG = pl.read_excel(Path("data/external/biolog_to_bigg_dict.xlsx"))

if __name__=="__main__":
    # Read data
    dark_data_path = Path("data/3_results/biolog/dark/growth_decision.xlsx")
    light_data_path = Path("data/3_results/biolog/light/growth_decision.xlsx")
    dark_data = pl.read_excel(
        dark_data_path,
        sheet_name="decision"
    )
    light_data = pl.read_excel(
        light_data_path,
        sheet_name="decision"
    )

    light_data = light_data.with_columns(condition=pl.lit("light"))
    dark_data = dark_data.with_columns(condition=pl.lit("dark"))

    data_df = pl.concat([light_data, dark_data])
    # Add bigg identifiers to the data
    data_df = data_df.join(
        BIOLOG_TO_BIGG.select("plate", "metabolite", "peptide", "new_bigg", "bigg", "metacyc"),
        on=["plate", "metabolite"]
    )
