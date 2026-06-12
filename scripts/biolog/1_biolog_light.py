"""Analisis of biolog data."""
from pathlib import Path

import numpy as np
import polars as pl
from labutils.models.growth import Gompertz
from labutils.tecan import load_tecan_scan_dir
from scipy.stats import ttest_ind
from xlsxwriter import Workbook

if __name__=="__main__":
    # Load the data
    data_dir = Path("data/0_raw/E260123/optical_density")
    files = [f for f in data_dir.iterdir() if f.is_file()]
    annot_cols = ["key", "condition"]
    data_df = load_tecan_scan_dir(files, annot_cols)
    data_df = data_df.with_columns(pl.col("sheet").str.slice(3, 2))
    data_df = data_df.rename({"key": "plate", "sheet": "replicate"})

    # Correct the absorbance with the pathlenght
    gpb = data_df.group_by("plate", "replicate", "well", "date")
    df_to_concat = []
    for name, data in gpb:
        abs980 = data.filter(pl.col("wavelength")==980)["absorbance"][0]
        abs900 = data.filter(pl.col("wavelength")==900)["absorbance"][0]
        volume = (abs980 - abs900) / 0.0010018
        pathlenght = (0.0055 * volume) + 0.0128
        corrected = (
            data
            .with_columns(
                (pl.col("absorbance") / pathlenght)
                .alias("abs_corrected"),
            )
        )
        df_to_concat.append(corrected)

    data_df = pl.concat(df_to_concat).sort("date", "plate", "well")

    # Fit Gompertz model
    x0 = [0.5, 0.1, 1]
    xmin = [0.01, 0.001, 0.01]
    model = Gompertz(x0, xmin)
    time = "days"

    query = data_df.filter(pl.col("wavelength")==680)
    gpb_cols = ["plate", "replicate", "well"]
    gpb = query.group_by(gpb_cols)

    data_list = []
    results_list = []
    for name, data in gpb:
        data = data.sort("date")
        t = data[time].to_numpy()
        # Log transform the absorbance
        y_min = data["abs_corrected"].min()
        y = data["abs_corrected"].to_numpy()
        y_log = np.log(y / y_min)

        # Fit
        fit = model.fit(t, y_log)
        error = fit.eval_uncertainty(sigma=0.9545, t=t)

        # Append best fit and log absorbance
        fitted_data = (
            data
            .with_columns(
                pl.lit(pl.Series(y_log)).alias("y_log"),
                pl.lit(pl.Series(fit.best_fit)).alias("y_pred"),
                pl.lit(pl.Series(error).alias("error")),
            )
        )
        data_list.append(fitted_data)

        new_dict = {
            "K": fit.params["K"].value,
            "mu": fit.params["mu"].value,
            "L": fit.params["L"].value,
            "K_err": fit.params["K"].stderr,
            "mu_err": fit.params["mu"].stderr,
            "L_err": fit.params["L"].stderr,
            "R2": fit.rsquared,
        }
        new_dict.update(dict(zip(gpb_cols, name)))
        results_list.append(new_dict)

    data_df = pl.concat(data_list)
    results_df = pl.DataFrame(results_list, orient="row")

    # Growth decision based on growth rates
    results_list = []
    gpb_cols = ["plate", "replicate"]
    gpb = results_df.group_by(gpb_cols)
    for name, data in gpb:
        # Eval if mu_well is larger that mu_ref + 1 * stddev
        mu_ref = data.filter(pl.col("well")=="A1")["mu"][0]
        std_ref = data.filter(pl.col("well")=="A1")["mu_err"][0]
        result = (
            data
            .with_columns(
                (pl.col("mu") > (mu_ref + std_ref)).alias("growth")
            )
        )
        results_list.append(result)
    results_df = pl.concat(results_list)

    # Growth is positive if True in more than 2 replicates
    cols = ["plate", "well"]
    gpb = results_df.group_by(cols)
    results_list = []
    for name, data in gpb:
        growth = data["growth"].sum() >= 2
        new_dict = {"growth": growth, "method": "mu"}
        new_dict.update(dict(zip(cols, name)))
        results_list.append(new_dict)

    decision_growth = pl.DataFrame(results_list, orient="row")

    # Growth decision based on p-vals
    # Growth decision based on p-values and log2fc
    p_thr = 0.05
    log2_fc_thr = 0.5

    gpb = results_df.group_by("plate")
    df_to_concat = []
    for name, data in gpb:
        # Get mu vals for well A1
        mu_ref = data.filter(pl.col("well")=="A1")["mu"].to_numpy()
        gpb_2 = data.group_by("well")
        for well, data_2 in gpb_2:
            # Get mu values for every other well
            mu_test = data_2["mu"].to_numpy()
            # Welch's t-test and get log2(test/ref)
            stat, p = ttest_ind(mu_ref, mu_test, equal_var=False)
            log2_fc = np.log2(np.mean(mu_test) / np.mean(mu_ref))
            # Save the data
            new_dict = {
                "plate": name[0],
                "well": well,
                "p-val": p,
                "log10_p": -np.log10(p),
                "log2_fc": log2_fc,
                # growth is True if p-val and log2FC are above the thresholds
                "method": "pvalue",
                "growth": (p < p_thr) and (log2_fc > log2_fc_thr),
            }
            df_to_concat.append(pl.DataFrame(new_dict))
    decision_pval = pl.concat(df_to_concat)

    # True in both methods
    a = decision_growth.filter(pl.col("growth")).select("plate", "well")
    b = decision_pval.filter(pl.col("growth")).select("plate", "well")
    decision_both = (
        b
        .join(a, on=["plate", "well"], how="inner")
        .with_columns(pl.lit(True).alias("growth"))
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
        .group_by("plate", "well")
        .agg(
            pl.col("mu").mean(),
            (pl.col("mu").std() + pl.col("mu_err").pow(2).sum() / 9).sqrt().alias("mu_err")
        )
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
        results_df.write_excel(
            wb,
            worksheet="fit",
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
