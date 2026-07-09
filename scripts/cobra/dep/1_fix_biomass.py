"""Normileze bimass reaction."""

from pathlib import Path

import cobra
import polars as pl
from memote.suite.tests.test_biomass import test_biomass_consistency

def fix_biomass_consistency(model, biomass_rxn_id):
    print(f"Testing reaction: {biomass_rxn_id}")
    biomass_rxn = model.reactions.get_by_id(biomass_rxn_id)
    mets_in_biomass = biomass_rxn.metabolites
    try:
        test_biomass_consistency(model, biomass_rxn_id)
    except AssertionError as e:
        print(e)
        mets_info = [(met.id, met.formula_weight, coeff) for met, coeff in mets_in_biomass.items()]
        mets_info_df = pl.DataFrame(
            mets_info,
            orient="row",
            schema=["met_id", "formula_weight", "coeff"]
        )
        new_coeffs = normalize_biomass_rxn(mets_info_df)
        biomass_rxn.subtract_metabolites(mets_in_biomass)
        biomass_rxn.add_metabolites(new_coeffs)
    else:
        print("Reaction is ok")

def normalize_biomass_rxn(mets_info_df):
    molecular_mass = 0
    for row in mets_info_df.iter_rows(named=True):
        molecular_mass += row["formula_weight"] * -row["coeff"]

    new_coeffs = {}
    new_mass = 0
    for row in mets_info_df.iter_rows(named=True):
        coeff = 1000 * row["coeff"] / molecular_mass
        new_mass += row["formula_weight"] * -coeff
        new_coeffs[row["met_id"]] = coeff

    print(f"The current molar mass is {molecular_mass / 1000}  mmol / g[CDW] / h")
    print(f"The new molar mass is {new_mass / 1000}  mmol / g[CDW] / h")

    return new_coeffs

if __name__=="__main__":
    model_dir = Path("models/draft/v0.0.1")
    model_file = model_dir / "hlacustris.xml"
    model = cobra.io.read_sbml_model(model_file)

    print("These mets don't have formula weight:\n")
    for met in model.reactions.BIOMASS_hlacus_auto.metabolites:
        if not met.formula:
            print(met)

    new_model = model.copy()
    biomass_rxn_id = "BIOMASS_hlacus_auto"
    fix_biomass_consistency(new_model, biomass_rxn_id)

    sol_old = model.slim_optimize()
    sol_new = new_model.slim_opitmize()
    print(f"Growth rate old: {sol_old}")
    print(f"Growth rate new: {sol_new}")

    cobra.io.write_sbml_model(new_model, "models/draft/v0.0.2")
