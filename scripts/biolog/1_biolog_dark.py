"""Analisis of biolog data."""
from pathlib import Path
from zipfile import ZipFile

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

def load_biolog_data(input_dir, plates, replicates):

    dt_expr = pl.col("date") - pl.col("date").first()
    ts_expr = dt_expr.dt.total_seconds()
    data_list = []
    for p in plates:
        for r in replicates:
            out_dir = input_dir / f"{p}/{r}"
            file = out_dir / "PM_Curation_Raw_HAL.zip"
            if not file.exists():
                print(f"FileNotFoundError: {str(file)}")
                continue
            raw_data_dir = out_dir / "HAL_RawReads"
            with ZipFile(file, "r") as zip:
                zip.extractall(out_dir)

            files = [f for f in raw_data_dir.iterdir() if f.is_file()]
            for raw_file in files:
                plate_df = pl.read_csv(raw_file)
                plate_df = (
                    plate_df
                    .with_columns(
                        pl.lit(p).alias("plate"),
                        pl.lit(r).alias("replicate"),
                        pl.col("Read At")
                        .str.to_datetime("%Y-%m-%dT%H:%M:%S%.fZ")
                        .alias("date"),
                    )
                    .sort("date")
                    .with_columns(
                        ts_expr.alias("seconds"),
                        (ts_expr / 60).alias("minutes"),
                        (ts_expr / 3600).alias("hours"),
                        (ts_expr / (24 * 3600)).alias("days"),
                    )
                    .drop("Read At", "Id", "PlateId", "Actual Temperature Celsius", "Target Temperature Celsius")
                )
                data_list.append(plate_df)

    data_df = pl.concat(data_list).rename({"Wavelength": "wavelength"})
    return (
        data_df
        .unpivot(
            index=["plate", "replicate", "wavelength", "date", "seconds", "minutes", "hours", "days"],
            variable_name="well",
            value_name="absorbance",
        )
    )

def down_sample_df(raw_df, samples_pers_day):
    h_to_keep = [int(h) for h in np.linspace(0, 24, samples_pers_day)]
    gpb_cols = gpb_cols = ["plate", "replicate", "well", "wavelength", "ordinal_day"]
    gpb = (
        raw_df
        .with_columns(
            pl.col("date").dt.ordinal_day().alias("ordinal_day"),
            pl.col("date").dt.hour().alias("ordinal_hour")
        )
        .group_by(gpb_cols, maintain_order=True)
    )

    dfs_to_append = []
    for name, data in gpb:
        new_df = (
            data
            .filter(pl.col("ordinal_hour").is_in(h_to_keep))
            .unique("ordinal_hour", keep="first")
        )
        dfs_to_append.append(new_df)

    return pl.concat(dfs_to_append)

if __name__=="__main__":
    # Load the data
    input_dir = "data/0_raw/E260123_biolog_dark"
    output_dir = "data/2_processed/E260123_biolog_dark"
    directory = {
        "input": Path(input_dir),
        "output": Path(output_dir),
    }
    directory["output"].mkdir(parents=True, exist_ok=True)

    plates = ["PM1", "PM2", "PM3"]
    replicates = ["R1", "R2", "R3"]
    raw_df = load_biolog_data(directory["input"], plates, replicates)

    # down sample the data to just N evenly spaced samples a day
    data_df = down_sample_df(raw_df, 2)

    # Fit Gompertz model
    x0 = [0.5, 0.1, 1]
    xmin = [0.01, 0.001, 0.0001]
    time = "days"
    gpb_cols = ["plate", "well"]
    control = "A01"
    data_df, results_df, fits = fit_groups_to_gompertz(
        data_df, 590, x0, xmin, gpb_cols, ["well", control], "days", "absorbance"
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
        mu_ref = data.filter(pl.col("well")==control)["mu"][0]
        std_ref = data.filter(pl.col("well")==control)["mu_err"][0]
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
        mu_ref = data.filter(pl.col("well")==control)["mu"][0]
        gpb_2 = data.group_by("well")
        for well, data_2 in gpb_2:
            # Get F-test and p-val
            fit_ref, _ = fits[f"{name[0]}_{control}"]
            fit_e, fit_r = fits[f"{name[0]}_{well[0]}"]
            chi2_r = fit_r.chisqr
            chi2_e = fit_e.chisqr + fit_ref.chisqr
            DF_r = len(fit_r.best_fit) - len(fit_r.params)
            DF_e = 2 * len(fit_e.best_fit) - 2 * len(fit_e.params)
            f_val = ((chi2_r - chi2_e) / (DF_r - DF_e)) / (chi2_e / DF_e)
            p_val = f.sf(f_val, (DF_r - DF_e), DF_e)

            # Get log2FC
            mu_test = data_2["mu"][0]
            log2_fc = np.log2(mu_test / mu_ref) if mu_test > 0 else 0.0
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
    biolog_map = (
        biolog_map
        .select("plate", "Well", "metabolite")
        .rename({"Well": "well"})
    )

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

    out_file = directory["output"] / "growth_decision.xlsx"
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
