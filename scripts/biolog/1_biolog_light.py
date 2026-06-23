"""Analisis of biolog data."""
from pathlib import Path

import numpy as np
import polars as pl
from labutils.models.growth import (
    Gompertz,
    fit_groups_to_gompertz,
)
from labutils.tecan import (
    load_tecan_scan_dir,
    correct_by_pathlenght,
)
from scipy.stats import ttest_ind_from_stats, f
from xlsxwriter import Workbook


if __name__=="__main__":
    # Load the data
    data_dir = Path("data/0_raw/E260123_biolog_light/optical_density")
    files = [f for f in data_dir.iterdir() if f.is_file()]
    annot_cols = ["key", "condition"]
    data_df = load_tecan_scan_dir(files, annot_cols)
    data_df = data_df.with_columns(pl.col("sheet").str.slice(3, 2))
    data_df = data_df.rename({"key": "plate", "sheet": "replicate"})

    # Correct the absorbance with the pathlenght
    gpby_cols = ["plate", "replicate", "well", "date"]
    data_df = correct_by_pathlenght(data_df, gpby_cols)



    # Fit Gompertz model
    x0 = [0.5, 0.1, 1]
    xmin = [0.01, 0.001, 0.0001]
    time = "days"
    gpb_cols = ["plate", "well"]
    data_df, results_df, fits = fit_groups_to_gompertz(
        data_df, 680, x0, xmin, gpb_cols, ["well", "A1"], "days"
    )

    # Growth decision based on growth rates
    decision_growth = (
        results_df
        .group_by(["plate", "well"])
        .agg(
            pl.col("mu").mean(),
            (pl.col("mu").std() + pl.col("mu_err").pow(2).sum() / 9)
            .sqrt().alias("mu_err"),
        )
    )

    gpb = results_df.group_by(["plate"])
    results_list = []
    for name, data in gpb:
        # Eval if mu_well is larger that mu_ref + 1 * stddev
        mu_ref = data.filter(pl.col("well")=="A1")["mu"][0]
        std_ref = data.filter(pl.col("well")=="A1")["mu_err"][0]
        result = (
            data
            .with_columns(
                ((pl.col("mu") - pl.col("mu_err")) > (mu_ref)).alias("growth")
            )
        )
        results_list.append(result)
    decision_growth = pl.concat(results_list)

    # Growth decision based on p-vals
    # Growth decision based on p-values and log2fc
    p_thr = 0.05
    log2_fc_thr = 1.5
    gpb = results_df.group_by("plate")
    df_to_concat = []
    for name, data in gpb:
        # Get mu vals for well A1
        mu_ref = data.filter(pl.col("well")=="A1")["mu"][0]
        gpb_2 = data.group_by("well")
        for well, data_2 in gpb_2:
            # Get F-test and p-val
            fit_ref, _ = fits[f"{name[0]}_A1"]
            fit_e, fit_r = fits[f"{name[0]}_{well[0]}"]
            chi2_r = fit_r.chisqr
            chi2_e = fit_e.chisqr + fit_ref.chisqr
            p = len(fit_e.params)
            N = 2*len(fit_e.best_fit)
            f_val = ((chi2_r - chi2_e) / p) / (chi2_e / (N - 2*p))
            p_val = f.sf(f_val, p, N - 2*p)

            # Get log2FC
            mu_test = data_2["mu"][0]
            log2_fc = np.log2(mu_test / mu_ref)
            # Save the data
            new_dict = {
                "plate": name[0],
                "well": well[0],
                "p_val": p_val,
                "f_val": f_val,
                "log10_p": -np.log10(p_val),
                "log2_fc": log2_fc,
                # growth is True if p-val and log2FC are above the thresholds
                "growth": (p_val < p_thr) and (log2_fc > log2_fc_thr),
            }
            df_to_concat.append(pl.DataFrame(new_dict))
    decision_pval = pl.concat(df_to_concat)

    # True in both methods
    a = decision_growth.filter(pl.col("growth")).select("plate", "well").sort("plate", "well")
    b = decision_pval.filter(pl.col("growth")).select("plate", "well")
    decision_both = (
        b
        .join(a, on=["plate", "well"], how="inner")
        .with_columns(pl.lit(True).alias("growth"))
        .sort("plate", "well")
    )

    biolog_map = pl.read_excel("data/external/biolog_map.xlsx")
    biolog_map = biolog_map.select("plate", "well", "metabolite")

    # Save results
    data_df = (
        data_df
        .join(biolog_map, on=["plate", "well"], how="left")
        .sort("plate", "well", "date")
    )
    results_df = (
        results_df
        .join(biolog_map, on=["plate", "well"], how="left")
        .sort("plate", "well")
    )
    decision_growth = (
        decision_growth
        .join(biolog_map, on=["plate", "well"], how="left")
        .sort("plate", "well")
    )
    decision_pval = (
        decision_pval
        .join(biolog_map, on=["plate", "well"], how="left")
        .sort("plate", "well")
    )
    decision_both = (
        decision_both
        .join(biolog_map, on=["plate", "well"], how="left")
        .sort("plate", "well")
    )
    missing = (
        biolog_map
        .join(decision_both, on=["plate", "well"], how="anti")
        .with_columns(pl.lit(False).alias("growth"))
        .select("plate", "well", "growth", "metabolite")
    )
    decision_final = pl.concat([decision_both, missing])
    query = (
        results_df
        .select("plate", "well", "mu", "mu_err")
        .join(decision_final, on=["plate", "well"])
    )

    out_dir = Path("data/2_processed/E260123_biolog_light/")
    out_dir.mkdir(parents=True, exist_ok=True)
    out_file = out_dir / "growth_decision.xlsx"
    with Workbook(out_file) as wb:
        data_df.write_excel(
            wb,
            worksheet="data",
        )
        decision_growth.write_excel(
            wb,
            worksheet="growth",
        )
        decision_pval.write_excel(
            wb,
            worksheet="pval",
        )
        query.write_excel(
            wb,
            worksheet="all",
        )
