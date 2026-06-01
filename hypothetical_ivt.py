import xarray as xr
import numpy as np
import os
import gc

# =========================
# Settings
# =========================
g = 9.81

# IVT layer in Pa
# Use 30000-100000 for 300-1000 hPa
p_top = 30000.0
p_bot = 100000.0

years = range(1979, 2026)

q_dir = "/pscratch/sd/s/slubis/ERA5_Data/SH/daily_mean"
u_dir = "/pscratch/sd/s/slubis/ERA5_Data/U/daily_mean"
v_dir = "/pscratch/sd/s/slubis/ERA5_Data/V/daily_mean"

q_clim_file = "./daily_climatology/q_daily_clim.nc"
u_clim_file = "./daily_climatology/u_daily_clim.nc"
v_clim_file = "./daily_climatology/v_daily_clim.nc"

out_dir = "./IVT_fixed_variable"
os.makedirs(out_dir, exist_ok=True)


# =========================
# Helper functions
# =========================
def get_pressure_name(ds):
    for name in ["plev", "level", "lev", "pressure"]:
        if name in ds.coords:
            return name
    raise ValueError("Cannot find pressure coordinate name.")


def pressure_layer_thickness(p, p_top, p_bot):
    """
    Approximate pressure-layer thickness dp for pressure levels.
    Pressure must be in Pa.
    """
    p = xr.DataArray(p)

    reverse_back = False
    if p[0] > p[-1]:
        p = p[::-1]
        reverse_back = True

    p_vals = p.values

    bounds = np.zeros(len(p_vals) + 1)
    bounds[1:-1] = 0.5 * (p_vals[:-1] + p_vals[1:])
    bounds[0] = p_top
    bounds[-1] = p_bot

    dp_vals = np.diff(bounds)

    dp = xr.DataArray(
        dp_vals,
        dims=p.dims,
        coords=p.coords,
        name="dp"
    )

    if reverse_back:
        dp = dp[::-1]

    return dp


def vertical_integral(field, dp, plev_name):
    return (field * dp / g).sum(dim=plev_name)


def magnitude(x, y):
    return np.sqrt(x**2 + y**2)


def compute_ivt(q, u, v, dp, plev_name):
    """
    Compute IVT magnitude from q, u, and v.
    """
    IVTx = vertical_integral(q * u, dp, plev_name)
    IVTy = vertical_integral(q * v, dp, plev_name)
    IVT = magnitude(IVTx, IVTy)
    return IVT


def match_daily_climatology(data, clim):
    """
    Match daily climatology to actual data time.

    Works for climatology with 365 or 366 time steps.
    Avoids xarray coordinate conflict by using numpy indices.
    """

    nclim = clim.sizes["time"]

    if nclim not in [365, 366]:
        raise ValueError(f"Climatology must have 365 or 366 days, found {nclim}")

    # Climatology has leap day
    if nclim == 366:
        doy_index = (data["time"].dt.dayofyear - 1).values
        clim_daily = clim.isel(time=doy_index)
        clim_daily = clim_daily.assign_coords(time=data["time"])
        return data, clim_daily

    # Climatology has no leap day
    is_feb29 = (data["time"].dt.month == 2) & (data["time"].dt.day == 29)

    if bool(is_feb29.any()):
        data = data.sel(time=~is_feb29)

    doy_index = (data["time"].dt.dayofyear - 1).values

    after_feb28_in_leap_year = (
        data["time"].dt.is_leap_year
        & (data["time"].dt.dayofyear > 59)
    )

    doy_index = doy_index - after_feb28_in_leap_year.values.astype(int)

    clim_daily = clim.isel(time=doy_index)
    clim_daily = clim_daily.assign_coords(time=data["time"])

    return data, clim_daily


# =========================
# Open climatology lazily
# =========================
chunks = {
    "time": 30,
    "plev": -1,
    "lat": 45,
    "lon": 90,
}

qclim_ds = xr.open_dataset(q_clim_file, chunks=chunks)
uclim_ds = xr.open_dataset(u_clim_file, chunks=chunks)
vclim_ds = xr.open_dataset(v_clim_file, chunks=chunks)

qclim = qclim_ds["var133"]
uclim = uclim_ds["var131"]
vclim = vclim_ds["var132"]

plev_name = get_pressure_name(uclim_ds)

# Pressure coordinate
p = uclim[plev_name]

# Convert hPa to Pa if needed
if p.max() < 2000:
    p_pa = p * 100.0
else:
    p_pa = p

# Select pressure range
plev_sel = p[(p_pa >= p_top) & (p_pa <= p_bot)]

qclim = qclim.sel({plev_name: plev_sel})
uclim = uclim.sel({plev_name: plev_sel})
vclim = vclim.sel({plev_name: plev_sel})

p = uclim[plev_name]

if p.max() < 2000:
    p = p * 100.0

dp = pressure_layer_thickness(p, p_top=p_top, p_bot=p_bot)
dp = dp.rename({dp.dims[0]: plev_name})

print("Using pressure levels:")
print(uclim[plev_name].values)

print("Climatology time lengths:")
print("qclim:", qclim.sizes["time"])
print("uclim:", uclim.sizes["time"])
print("vclim:", vclim.sizes["time"])


# =========================
# Loop over years
# =========================
for year in years:
    print(f"Processing year: {year}")

    q_file = f"{q_dir}/q.{year}.nc"
    u_file = f"{u_dir}/u.{year}.nc"
    v_file = f"{v_dir}/v.{year}.nc"

    q_ds = xr.open_dataset(q_file, chunks=chunks)
    u_ds = xr.open_dataset(u_file, chunks=chunks)
    v_ds = xr.open_dataset(v_file, chunks=chunks)

    q = q_ds["var133"].sel({plev_name: plev_sel})
    u = u_ds["var131"].sel({plev_name: plev_sel})
    v = v_ds["var132"].sel({plev_name: plev_sel})

    # =========================
    # Match daily climatology
    # =========================
    q, q_clim_daily = match_daily_climatology(q, qclim)
    u, u_clim_daily = match_daily_climatology(u, uclim)
    v, v_clim_daily = match_daily_climatology(v, vclim)

    # Ensure q, u, and v have the same time axis
    common_time = np.intersect1d(
        np.intersect1d(q["time"].values, u["time"].values),
        v["time"].values
    )

    q = q.sel(time=common_time)
    u = u.sel(time=common_time)
    v = v.sel(time=common_time)

    q_clim_daily = q_clim_daily.sel(time=common_time)
    u_clim_daily = u_clim_daily.sel(time=common_time)
    v_clim_daily = v_clim_daily.sel(time=common_time)

    print(f"Number of daily time steps used for {year}: {q.sizes['time']}")

    # =========================
    # Fixed-variable IVT calculations
    # =========================

    # Total IVT:
    # actual moisture and actual winds
    IVT = compute_ivt(q, u, v, dp, plev_name)

    # Dynamic IVT:
    # moisture held fixed to daily climatology, winds vary
    IVT_dyn = compute_ivt(q_clim_daily, u, v, dp, plev_name)

    # Thermodynamic IVT:
    # winds held fixed to daily climatology, moisture varies
    IVT_thermo = compute_ivt(q, u_clim_daily, v_clim_daily, dp, plev_name)

    # =========================
    # Output only requested variables
    # =========================
    ds_out = xr.Dataset(
        {
            "IVT": IVT.astype("float32"),
            "IVT_dyn": IVT_dyn.astype("float32"),
            "IVT_thermo": IVT_thermo.astype("float32"),
        }
    )

    ds_out["IVT"].attrs["long_name"] = (
        "Integrated vapor transport using actual moisture and actual winds"
    )

    ds_out["IVT_dyn"].attrs["long_name"] = (
        "Fixed-moisture IVT using daily climatological moisture and actual winds"
    )

    ds_out["IVT_thermo"].attrs["long_name"] = (
        "Fixed-wind IVT using actual moisture and daily climatological winds"
    )

    for var in ds_out.data_vars:
        ds_out[var].attrs["units"] = "kg m-1 s-1"

    ds_out.attrs["description"] = (
        "Fixed-variable IVT sensitivity calculations. "
        "IVT is computed from actual q, u, and v. "
        "IVT_dyn is computed with q held fixed to its daily climatological value. "
        "IVT_thermo is computed with u and v held fixed to their daily climatological values. "
        "Anomalies should be calculated separately by subtracting each variable's own daily climatology."
    )

    encoding = {
        var: {
            "zlib": True,
            "complevel": 4,
            "dtype": "float32",
            "_FillValue": np.float32(np.nan),
        }
        for var in ds_out.data_vars
    }

    out_file = f"{out_dir}/IVT_fixed_variable_{year}.nc"

    print(f"Writing {out_file}")
    ds_out.to_netcdf(out_file, encoding=encoding)

    print(f"Saved: {out_file}")

    q_ds.close()
    u_ds.close()
    v_ds.close()

    del q_ds, u_ds, v_ds
    del q, u, v
    del q_clim_daily, u_clim_daily, v_clim_daily
    del IVT, IVT_dyn, IVT_thermo
    del ds_out

    gc.collect()


# =========================
# Close climatology files
# =========================
qclim_ds.close()
uclim_ds.close()
vclim_ds.close()
