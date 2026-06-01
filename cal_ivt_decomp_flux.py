import xarray as xr
import numpy as np
import os

# =========================
# Settings
# =========================
g = 9.81
p_top = 100.0       # Pa
p_bot = 100000.0    # Pa

years = range(1979, 1981)

q_dir = "/pscratch/sd/s/slubis/ERA5_Data/SH/daily_mean"
u_dir = "/pscratch/sd/s/slubis/ERA5_Data/U/daily_mean"
v_dir = "/pscratch/sd/s/slubis/ERA5_Data/V/daily_mean"

# Daily climatology files containing q, u, v daily climatology
q_clim_file = "./daily_climatology/q_daily_clim.nc"
u_clim_file = "./daily_climatology/u_daily_clim.nc"
v_clim_file = "./daily_climatology/v_daily_clim.nc"

out_dir = "./IVT_decomposition"
os.makedirs(out_dir, exist_ok=True)


# =========================
# Helper functions
# =========================
def get_pressure_name(ds):
    for name in ["plev", "level", "lev", "pressure"]:
        if name in ds.coords:
            return name
    raise ValueError("Cannot find pressure coordinate name.")


def pressure_layer_thickness(p, p_top=100.0, p_bot=100000.0):
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


def match_daily_climatology(data, clim):
    """
    Match daily climatology to the actual data time axis.

    Handles:
    1. climatology with 366 days
    2. climatology with 365 days
    3. leap years in the input data
    4. xarray time-coordinate conflict by using numpy index arrays
    """

    nclim = clim.sizes["time"]

    if nclim not in [365, 366]:
        raise ValueError(
            f"Climatology must have 365 or 366 time steps, but found {nclim}"
        )

    # If climatology has 366 days, use normal day-of-year indexing
    if nclim == 366:
        doy_index = (data["time"].dt.dayofyear - 1).values

        clim_daily = clim.isel(time=doy_index)
        clim_daily = clim_daily.assign_coords(time=data["time"])

        return data, clim_daily

    # If climatology has 365 days, remove Feb 29 from data
    # and shift leap-year dates after Feb 28 backward by one index.
    if nclim == 365:
        is_feb29 = (data["time"].dt.month == 2) & (data["time"].dt.day == 29)

        if bool(is_feb29.any()):
            data = data.sel(time=~is_feb29)

        doy_index = (data["time"].dt.dayofyear - 1).values

        # For leap years, after Feb 28, subtract one day from index
        after_feb28_in_leap_year = (
            data["time"].dt.is_leap_year
            & (data["time"].dt.dayofyear > 59)
        )

        doy_index = doy_index - after_feb28_in_leap_year.values.astype(int)

        clim_daily = clim.isel(time=doy_index)
        clim_daily = clim_daily.assign_coords(time=data["time"])

        return data, clim_daily


# =========================
# Read daily climatology
# =========================
qclim_ds = xr.open_dataset(q_clim_file)
uclim_ds = xr.open_dataset(u_clim_file)
vclim_ds = xr.open_dataset(v_clim_file)

qclim = qclim_ds["var133"]
uclim = uclim_ds["var131"]
vclim = vclim_ds["var132"]

print("Climatology time lengths:")
print("qclim:", qclim.sizes["time"])
print("uclim:", uclim.sizes["time"])
print("vclim:", vclim.sizes["time"])

plev_name = get_pressure_name(uclim_ds)
p = uclim[plev_name]

# Convert hPa to Pa if needed
if p.max() < 2000:
    p = p * 100.0

dp = pressure_layer_thickness(p, p_top=p_top, p_bot=p_bot)
dp = dp.rename({dp.dims[0]: plev_name})


# =========================
# Loop over years
# =========================
for year in years:
    print(f"Processing year: {year}")

    q_file = f"{q_dir}/q.{year}.nc"
    u_file = f"{u_dir}/u.{year}.nc"
    v_file = f"{v_dir}/v.{year}.nc"

    q_ds = xr.open_dataset(q_file)
    u_ds = xr.open_dataset(u_file)
    v_ds = xr.open_dataset(v_file)

    q = q_ds["var133"]
    u = u_ds["var131"]
    v = v_ds["var132"]

    plev_name = get_pressure_name(u_ds)

    # =========================
    # Match daily climatology by day of year
    # =========================
    q, q_clim_daily = match_daily_climatology(q, qclim)
    u, u_clim_daily = match_daily_climatology(u, uclim)
    v, v_clim_daily = match_daily_climatology(v, vclim)

    # Make sure q, u, and v have identical time axes after possible Feb 29 removal
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
    # Anomalies relative to daily climatology
    # =========================
    q_prime = q - q_clim_daily
    u_prime = u - u_clim_daily
    v_prime = v - v_clim_daily

    # =========================
    # Total IVT
    # =========================
    IVTx = vertical_integral(q * u, dp, plev_name)
    IVTy = vertical_integral(q * v, dp, plev_name)
    IVT = magnitude(IVTx, IVTy)

    # =========================
    # Vector decomposition components
    # =========================

    # Dynamic term: q_clim * V'
    IVTx_dyn = vertical_integral(q_clim_daily * u_prime, dp, plev_name)
    IVTy_dyn = vertical_integral(q_clim_daily * v_prime, dp, plev_name)

    # Thermodynamic/moisture term: q' * V_clim
    IVTx_thermo = vertical_integral(q_prime * u_clim_daily, dp, plev_name)
    IVTy_thermo = vertical_integral(q_prime * v_clim_daily, dp, plev_name)

    # Nonlinear term: q' * V'
    IVTx_nonlinear = vertical_integral(q_prime * u_prime, dp, plev_name)
    IVTy_nonlinear = vertical_integral(q_prime * v_prime, dp, plev_name)

    # Magnitudes of each vector component
    IVT_dyn = magnitude(IVTx_dyn, IVTy_dyn)
    IVT_thermo = magnitude(IVTx_thermo, IVTy_thermo)
    IVT_nonlinear = magnitude(IVTx_nonlinear, IVTy_nonlinear)

    # Optional: total anomaly vector reconstruction
    IVTx_clim = vertical_integral(q_clim_daily * u_clim_daily, dp, plev_name)
    IVTy_clim = vertical_integral(q_clim_daily * v_clim_daily, dp, plev_name)

    IVTx_anom = IVTx - IVTx_clim
    IVTy_anom = IVTy - IVTy_clim
    IVT_anom = magnitude(IVTx_anom, IVTy_anom)

    # This vector identity should hold:
    # IVTx_anom = IVTx_dyn + IVTx_thermo + IVTx_nonlinear
    # IVTy_anom = IVTy_dyn + IVTy_thermo + IVTy_nonlinear

    # =========================
    # Save variables
    # =========================
    ds_out = xr.Dataset(
        {
            "IVT": IVT,
            "IVT_anom": IVT_anom,

            "IVT_dyn": IVT_dyn,
            "IVT_thermo": IVT_thermo,
            "IVT_nonlinear": IVT_nonlinear,

            "IVTx": IVTx,
            "IVTy": IVTy,

            "IVTx_anom": IVTx_anom,
            "IVTy_anom": IVTy_anom,

            "IVTx_dyn": IVTx_dyn,
            "IVTy_dyn": IVTy_dyn,

            "IVTx_thermo": IVTx_thermo,
            "IVTy_thermo": IVTy_thermo,

            "IVTx_nonlinear": IVTx_nonlinear,
            "IVTy_nonlinear": IVTy_nonlinear,
        }
    )

    ds_out["IVT"].attrs["long_name"] = "Total integrated vapor transport magnitude"
    ds_out["IVT_anom"].attrs["long_name"] = "Magnitude of IVT vector anomaly"

    ds_out["IVT_dyn"].attrs["long_name"] = (
        "Magnitude of dynamic contribution to IVT anomaly: "
        "climatological moisture times anomalous wind"
    )
    ds_out["IVT_thermo"].attrs["long_name"] = (
        "Magnitude of thermodynamic contribution to IVT anomaly: "
        "anomalous moisture times climatological wind"
    )
    ds_out["IVT_nonlinear"].attrs["long_name"] = (
        "Magnitude of nonlinear contribution to IVT anomaly: "
        "anomalous moisture times anomalous wind"
    )

    ds_out["IVTx"].attrs["long_name"] = "Zonal component of total IVT"
    ds_out["IVTy"].attrs["long_name"] = "Meridional component of total IVT"

    ds_out["IVTx_anom"].attrs["long_name"] = "Zonal component of IVT anomaly"
    ds_out["IVTy_anom"].attrs["long_name"] = "Meridional component of IVT anomaly"

    ds_out["IVTx_dyn"].attrs["long_name"] = "Zonal dynamic IVT anomaly contribution"
    ds_out["IVTy_dyn"].attrs["long_name"] = "Meridional dynamic IVT anomaly contribution"

    ds_out["IVTx_thermo"].attrs["long_name"] = "Zonal thermodynamic IVT anomaly contribution"
    ds_out["IVTy_thermo"].attrs["long_name"] = "Meridional thermodynamic IVT anomaly contribution"

    ds_out["IVTx_nonlinear"].attrs["long_name"] = "Zonal nonlinear IVT anomaly contribution"
    ds_out["IVTy_nonlinear"].attrs["long_name"] = "Meridional nonlinear IVT anomaly contribution"

    for var in ds_out.data_vars:
        ds_out[var].attrs["units"] = "kg m-1 s-1"

    ds_out.attrs["description"] = (
        "Integrated vapor transport decomposition. "
        "The vector anomaly satisfies IVT_anom_vector = dynamic + thermodynamic + nonlinear. "
        "Magnitude terms are saved for convenience but should not be summed directly."
    )

    out_file = f"{out_dir}/IVT_decomposition_{year}.nc"
    ds_out.to_netcdf(out_file)

    print(f"Saved: {out_file}")

    q_ds.close()
    u_ds.close()
    v_ds.close()

qclim_ds.close()
uclim_ds.close()
vclim_ds.close()
