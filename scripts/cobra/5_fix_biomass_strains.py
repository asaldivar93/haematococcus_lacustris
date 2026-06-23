"""Fix non-growing model."""
from pathlib import Path

import cobra
import polars as pl
from labutils.cobra.biomass import fix_biomass_mw
from labutils.cobra.io import write_excel
from tqdm import tqdm


def find_difference(model_1, model_2):
    model_1_df = [(rxn.id, rxn.bounds, list(rxn.compartments))
        for rxn in model_1.reactions]
    model_1_df = pl.DataFrame(
        model_1_df,
        orient="row",
        schema=["rxn_id_r", "bounds_r", "compartments_r"],
    )

    model_2_df = [(rxn.id, rxn.bounds, list(rxn.compartments))
        for rxn in model.reactions]
    model_2_df = pl.DataFrame(
        model_2_df,
        orient="row",
        schema=["rxn_id_l", "bounds_l", "compartments_l"],
    )

    # Find reactions that are different in the reference and the model
    join_expr = ((pl.col("bounds_l") != pl.col("bounds_r"))
        | (pl.col("compartments_l") != pl.col("compartments_r")))
    different = (
        model_1_df
        .join_where(
            model_2_df,
            pl.col("rxn_id_r")==pl.col("rxn_id_l"),
            join_expr,
        )
    )

    return different

def find_essential_rxns(model, rxns_to_test):
    base_growth = model.slim_optimize()
    essential_rxns = []
    for rxn in tqdm(rxns_to_test):
        with model as m:
            r = m.reactions.get_by_id(rxn.id)
            m.remove_reactions([r])
            growth = m.slim_optimize()
        if abs(base_growth - growth) >= 1e-6:
            essential_rxns.append(rxn)
    return essential_rxns

if __name__=="__main__":
    # Reference model
    strain = "redball"
    model_file = Path(f"models/draft/v0.0.5/{strain}/{strain}.xml")
    reference = cobra.io.read_sbml_model(model_file)
    ref_rxns = {rxn.id for rxn in reference.reactions}

    # Load the model
    strain = "nies144"
    model_file = Path(f"models/draft/v0.0.5/{strain}/{strain}.xml")
    model = cobra.io.read_sbml_model(model_file)
    reactions_df = pl.read_excel(
        f"models/draft/v0.0.5/{strain}/{strain}.xlsx",
        sheet_name="reactions"
    )

    # Normalize biomass reaction
    biomass_rxns = ["BIOMASS_hlacus_auto", "BIOMASS_hlacus_mixo", "BIOMASS_hlacus_hetero"]
    for rxn_id in biomass_rxns:
        fix_biomass_mw(model, rxn_id)

    # Find reactions that are different in the reference and the model
    different = find_difference(reference, model)

    # Substitute reactions in the model with reactions from the reference
    rxns_to_remove = []
    rxns_to_add = []
    for rxn_id in different["rxn_id_l"]:
        ref = reference.reactions.get_by_id(rxn_id).copy()
        this = model.reactions.get_by_id(rxn_id)
        ref.gpr = this.gpr
        rxns_to_remove.append(this)
        rxns_to_add.append(ref)

    model.remove_reactions(rxns_to_remove)
    model.add_reactions(rxns_to_add)

    # Add all missing reactions model
    model_rxns = {rxn.id for rxn in model.reactions}
    missing_rxns = ref_rxns - model_rxns
    rxns_to_add = [reference.reactions.get_by_id(rxn_id).copy()
        for rxn_id in missing_rxns]
    model.add_reactions(rxns_to_add)

    model.slim_optimize()



    # Save the model
    out_dir = Path("models/draft/v0.0.6/nies144/")
    out_dir.mkdir(parents=True, exist_ok=True)
    out_file = out_dir / "nies144.xlsx"
    write_excel(model, out_file)

    out_file = out_dir / "nies144.xml"
    cobra.io.write_sbml_model(model, out_file)

    out_file = out_dir / "nies144.json"
    cobra.io.save_json_model(model, out_file)
