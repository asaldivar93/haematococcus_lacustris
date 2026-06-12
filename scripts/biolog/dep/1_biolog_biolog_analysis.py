import importlib.resources
import json
import os
from pathlib import Path


import numpy as np
import plotly.express as px
import plotly.graph_objects as go
import plotly.io as pio
import polars as pl
from pynumdiff.optimize import optimize
from pynumdiff.utils import evaluate
from tqdm import tqdm

from pynumdiff.finite_difference import finitediff
from pynumdiff.smooth_finite_difference import kerneldiff, butterdiff
from pynumdiff.polynomial_fit import splinediff, polydiff, savgoldiff
from pynumdiff.basis_fit import spectraldiff, rbfdiff
from pynumdiff.total_variation_regularization import tvrdiff, smooth_acceleration
from pynumdiff.kalman_smooth import rtsdiff, robustdiff
from pynumdiff.linear_model import lineardiff

methods_dict = {
    "finitediff": finitediff, "kerneldiff": kerneldiff, "splinediff": splinediff,
    "polydiff": polydiff, "savgoldiff": savgoldiff, "spectraldiff": spectraldiff,
    "rbfdiff": rbfdiff, "tvrdiff": tvrdiff, "smooth_acceleration": smooth_acceleration,
    "rtsdiff": rtsdiff, "robustdiff": robustdiff, "lineardiff": lineardiff,
}

pio.templates.default = "simple_white"

########## Load the data ##############
plates = ["PM1", "PM2", "PM3"]
replicates = ["R1", "R2", "R3"]

data_list = []
for p in plates:
    for r in replicates:
        out_dir = f"data/0_raw/biolog/E260123/{p}/{r}"
        file = f"{out_dir}/PM_Curation_Raw_HAL.zip"
        cmd = f"unzip {file} -d {out_dir}"
        os.system(cmd)

        try:
            file = os.listdir(f"{out_dir}/HAL_RawReads")[0]
        except FileNotFoundError as e:
            print(e)
        else:
            raw_data_file = Path(f"{out_dir}/HAL_RawReads/{file}")
            print(raw_data_file)
            plate_df = pl.scan_csv(raw_data_file)
            plate_df = (
                plate_df
                .with_columns(
                    pl.lit(p).alias("Plate"),
                    pl.lit(r).alias("Replicate"),
                    pl.col("Read At")
                    .str.to_datetime("%Y-%m-%dT%H:%M:%S%.fZ")
                    .alias("Date"),
                )
                .drop("Read At", "Id", "PlateId", "Actual Temperature Celsius", "Target Temperature Celsius")
                .sort("Date")
                .with_columns(Time=pl.col("Date") - pl.col("Date").first())
                .with_columns(Time=pl.col("Time").dt.total_seconds() / (3600))
            )
            data_list.append(plate_df)

data_df = pl.concat(data_list)
data_df = (
    data_df
    .unpivot(
        index=["Date","Time","Wavelength","Plate","Replicate"],
        variable_name="Well",
        value_name="abs"
    )
)

######### Do some plots ###########

# Plot data
q = (
    data_df
    .group_by("Wavelength", "Plate", "Replicate", "Well")
    .agg([
        (pl.col("abs").last() - pl.col("abs").first())
    ])
    .sort("abs", descending=True)
)

q = (
    data_df
    .filter(
        pl.col("Plate")=="PM3",
        pl.col("Replicate")=="R2",
        pl.col("Well")=="H05",
        pl.col("Wavelength")==590,
    )
    .sort("Time")
)

# Plot Data
px.scatter(
    q.collect(),
    x="Time",
    y="abs"
)

######### Select the numerical differentiation method ###########
# Plot Power spectrum
avg_dt = (q.collect()["Time"][1:] - q.collect()["Time"][0:-1]).mean()
time = q.collect()["Time"]
x = q.collect()["abs"].to_numpy()
X = np.fft.fft(x)
energy = 20*np.log10(np.abs(X))
freqs = np.fft.fftfreq(len(X), avg_dt)
energy = energy[freqs >=0]
freqs = freqs[freqs >= 0]

go.Figure().add_trace(
    go.Scatter(
        x=freqs, y=energy
    )
)

# from the plot choose a cutoff frequency


# Plot different numeric methods
cutoff_freq = 0.5
log_gamma = -1.6*np.log(cutoff_freq) -0.71*np.log(avg_dt) - 5.1
tvgamma = np.exp(log_gamma)

fig1 = go.Figure()
fig2 = go.Figure()
for name, method in tqdm(methods_dict.items()):
    params, val = optimize(method, x, avg_dt, tvgamma=tvgamma)
    x_hat, dxdt_hat = method(x, avg_dt, **params)
    fig1.add_trace(
        go.Scatter(x=time, y=x_hat, name=name)
    )
    fig2.add_trace(
        go.Scatter(x=time, y=dxdt_hat, name=name)
    )

fig2.show()
# From the dxdt_hat plot, choose a method to use

# Plot different cut_off freqs
freqs = [0.5, 0.75, 1, 1.5]
fig3 = go.Figure()
for f in freqs:
    log_gamma = -1.6*np.log(f) -0.71*np.log(avg_dt) - 5.1
    tvgamma = np.exp(log_gamma)
    method = methods_dict["splinediff"]
    params, val = optimize(method, x, avg_dt, tvgamma=tvgamma)
    x_hat, dxdt_hat = method(x, avg_dt, **params)
    fig3.add_trace(
        go.Scatter(x=time, y=dxdt_hat, name=f)
    )

######## Get instantaneus rates
######## This is the code to get the instantaneus rates #########
######## You can ignore most of the rest ######

cutoff_freq = 0.5
log_gamma = -1.6*np.log(cutoff_freq) -0.71*np.log(avg_dt) - 5.1
tvgamma = np.exp(log_gamma)

gpb = data_df.collect().group_by("Plate", "Replicate", "Well", "Wavelength")
method = methods_dict["splinediff"]
params, val = optimize(method, x, avg_dt, tvgamma=tvgamma)

def get_dxdt(group: pl.DataFrame) -> pl.DataFrame:
    x = group["abs"].to_numpy()
    x_hat, dxdt_hat = method(x, avg_dt, **params)
    return group.with_columns(pl.Series("dxdt_hat", dxdt_hat))

keys = ["Plate", "Replicate", "Well", "Wavelength"]

data_df = (
    data_df
    .group_by(keys)
    .map_groups(get_dxdt, schema=None)
    .collect()
    .sort(keys)
)

rates_df = (
    data_df
    .group_by(keys)
    .agg((pl.col("dxdt_hat").max() * 24).alias("mu (day-1)"))
    .sort(keys)
)

out_dir = Path("data/1_interim/biolog/E260123/")
out_dir.mkdir(parents=True, exist_ok=True)
data_df.write_csv(out_dir / "data_df.csv")
rates_df.write_csv(out_dir / "rates_df.csv")
