from pathlib import Path

import cobra
import polars as pl
from labutils.cobra.io import write_excel
from labutils.cobra.mfg import find_blocked_mets
from sklearn.metrics import balanced_accuracy_score
from tqdm import tqdm
from xlsxwriter import Workbook


def get_biologs_in_model(model, positives):
    mets_in_model = [met.id for met in model.metabolites]
    compartments = list(model.compartments)
    query = (
        positives
        .unique("metabolite")
        .select("metabolite", "bigg.metabolite", "peptide")
    )
    results = []
    for row in query.iter_rows(named=True):
        mets = row["bigg.metabolite"].split(",")
        for met_id in mets:
            comps = {c: f"{met_id}_{c}" for c in compartments
                if f"{met_id}_{c}" in mets_in_model}
            comps["bigg.metabolite"] = met_id
            comps["metabolite"] = row["metabolite"]
            results.append(comps)

    biologs_with_bigg = pl.DataFrame(results, orient="row")
    target_cols = [c for c in biologs_with_bigg.columns if c in compartments]
    query = (
        biologs_with_bigg
        .with_columns(
            pl.all_horizontal(pl.col(target_cols).is_null()).alias("all_null")
        )
        .select("metabolite", "bigg.metabolite", "all_null", *target_cols)
    )
    biologs_in_model = query.filter(~pl.col("all_null")).select(pl.exclude("all_null"))
    missing_biologs = query.filter(pl.col("all_null")).select(pl.exclude("all_null"))

    return biologs_in_model, missing_biologs

def simulate_row(model, plate, condition, met_id):
    m = model.copy()
    base_medium = m.medium
    pm_media = {
        "PM1": {k: v for k, v in base_medium.items() if k!="EX_co2_e"},
        "PM2": {k: v for k, v in base_medium.items() if k!="EX_co2_e"},
        "PM3": {k: v for k, v in base_medium.items() if k!="EX_no3_e"},
    }

    rxn_id = "EX_" + met_id
    medium = pm_media[plate].copy()
    medium[rxn_id] = 100
    match condition:
        case "dark":
            medium["EX_photonVis_e"]=0

    try:
        m.medium = medium
        has_exchange = True
    except KeyError:
        has_exchange = False
        met = m.metabolites.get_by_id(row["e"])
        m.add_boundary(met, type="exchange")
        m.medium = medium

    return m, m.optimize(), has_exchange

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
        keys_to_remove = ["date", "met_id", "formula", "name", "charge"]
        for key in keys_to_remove:
            met.pop(key)
        new_met.annotations = met

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
        keys_to_remove = ["date", "rxn_id", "name", "lower_bound",
            "upper_bound", "gpr", "subsystem", "gapfilling", "localized"]
        for key in keys_to_remove:
            row.pop(key)
        new_rxn.annotation = row

        rxns_to_add.append(new_rxn)
    return rxns_to_add

def update_boundary_rxns(model, boundary):
    for row in boundary.iter_rows(named=True):
        met = model.metabolites.get_by_id(row["met_id"])
        model.add_boundary(met, type=row["type"])

def update_model(model, updates_path):
    metabolites = pl.read_excel(updates_path, sheet_name="metabolites")
    #reactions = pl.read_excel(updates_path, sheet_name="reactions")
    #boundary = pl.read_excel(updates_path, sheet_name="boundary_rxns")

    mets_to_add = update_metabolites(metabolites)
    #rxns_to_add = update_reactions(reactions)
    #update_boundary_rxns(model, boundary)

    model.add_metabolites(mets_to_add)
    #model.add_reactions(rxns_to_add)


if __name__=="__main__":
    # Inputs
    model_path = Path("models/draft/v0.0.6/nies144/nies144.xml")
    worklist_path = Path("data/metrics.list")
    updates_path = Path("data/1_interim/curation/updates_to_model.xlsx")
    output_path = Path("data/2_processed/biolog_metrics.xlsx")

    # Load Model
    model = cobra.io.read_sbml_model(model_path)
    # Update the model
    update_model(model, updates_path)

    # Find and remove deadend and orphan metabolites
    rxns_to_remove = True
    removed_rxns = []
    removed_mets = []
    while rxns_to_remove:
        # orphans, deadends = find_blocked_mets(model)

        rxns_to_remove = []
        for met in model.metabolites:
            # met = model.metabolites.get_by_id(mid)
            if len(met.reactions) <= 1:
                rxns_to_remove.extend([*met.reactions])

        rxns_to_remove = set(rxns_to_remove)
        model.remove_reactions(rxns_to_remove)
        removed_rxns.extend(rxns_to_remove)

        mets_to_remove = [met for met in model.metabolites if len(met.reactions)==0]
        model.remove_metabolites(mets_to_remove)
        removed_mets.extend(mets_to_remove)

    genes_to_remove = [g.id for g in model.genes if not bool(g.reactions)]
    for gid in genes_to_remove:
        model.genes.remove(gid)
    write_excel(model, "models/draft/v0.0.7/nies144/nies144.xlsx")
    model.slim_optimize()
    model.metabolites.get_by_id("2oxoadp_m")
    # Load biolog mets info
    biolog_file = Path("data/external/biolog_map.xlsx")
    biolog_map = pl.read_excel(biolog_file)

    # Load biolog results
    worklist_df = pl.read_csv(worklist_path)
    results_list = []
    for row in worklist_df.iter_rows(named=True):
        new_data = pl.read_excel(Path(row["path"]), sheet_name="all")
        new_data = new_data.with_columns(pl.lit(row["condition"]).alias("condition"))
        results_list.append(new_data)
    biolog_data = pl.concat(results_list)
    query = (
        biolog_data
        .join(
            biolog_map.select("plate", "metabolite", "bigg.metabolite", "peptide"),
            on=["plate", "metabolite"]
        )
        .drop_nulls("bigg.metabolite")
    )

    # Get biolog mets already in the model
    biologs_in_model, missing_biologs = get_biologs_in_model(model, query)

    # Not in the model but positive
    query = (
        biolog_data
        .select("metabolite", "growth")
        .filter(pl.col("growth"))
    )
    missing_in_model = (
        missing_biologs
        .join(query, on="metabolite")
    )
    # In model but not in the external compartment
    missing_e = (
        biolog_data.select("metabolite")
        .join(biologs_in_model, on="metabolite")
        .filter(pl.col("e").is_null())
        .unique("metabolite")
    )
    to_simulate = (
        biolog_data
        .join(biologs_in_model, on="metabolite")
        .filter(~pl.col("e").is_null())
    )
    mssg_str = f"From {biolog_data.shape[0]} biologs \
    {to_simulate.shape[0]} are ready for simulation"
    print(mssg_str)

    # Do simulationsEX_no3_e
    tol = 1e-6
    sim_results = []
    it = to_simulate.iter_rows(named=True)
    for row in tqdm(it, total=to_simulate.shape[0]):
        m, sol, has_exchange = simulate_row(
            model, row["plate"], row["condition"], row["e"]
        )
        mu_pred = sol.objective_value
        y_pred = mu_pred > tol
        row["y_pred"] = y_pred
        row["mu_pred"] = mu_pred
        row["has_exchange"] = has_exchange
        row["solution"] = sol
        row["instance"] = m

        match row["y_pred"] == row["growth"]:
            case True:
                if row["growth"]:
                    row["type"] = "TP"
                else:
                    row["type"] = "TN"
            case False:
                if row["growth"]:
                    row["type"] = "FN"
                else:
                    row["type"] = "FP"

        sim_results.append(row)

    col_order = ["plate", "well", "condition", "metabolite", "bigg.metabolite",
                 "m", "h", "c", "x", "e", "mu", "mu_err", "mu_pred",
                 "growth", "y_pred", "type", "has_exchange", "solution", "instance"]

    results = pl.DataFrame(sim_results, orient="row").select(col_order)
    results = results.sort("plate",  "condition", "well").with_row_index("rid")

    query = (
        results
        .select(pl.exclude(["solution", "instance"]))
    )
    with Workbook(output_path) as wb:
        missing_in_model.write_excel(wb, worksheet="missing")
        missing_e.write_excel(wb, worksheet="missing_e")
        query.write_excel(wb, worksheet="results")

    balanced_accuracy_score(results["growth"], results["y_pred"])
    # Evaluate the solutions
    solutions = dict(zip(results["rid"], results["solution"]))
    instances= dict(zip(results["rid"], results["instance"]))
    sol = solutions[0]
    m = instances[0]
    m.summary(solution=sol)
    met_id="g6p_A_c"
    met = m.metabolites.get_by_id(met_id)
    met.summary(solution=sol)
