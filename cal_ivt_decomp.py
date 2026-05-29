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

q_dir = "/pscratch/sd/s/slubis/ERA5_Data/SH/daily"
u_dir = "/pscratch/sd/s/slubis/ERA5_Data/U/daily"
v_dir = "/pscratch/sd/s/slubis/ERA5_Data/V/daily"

# Daily climatology files containing q, u, v daily climatology
q_clim_file = "/path/to/q_daily_clim.nc"
u_clim_file = "/path/to/u_daily_clim.nc"
v_clim_file = "/path/to/v_daily_clim.nc"

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


# =========================
# Read daily climatology
# =========================
qclim_ds = xr.open_dataset(q_clim_file)
uclim_ds = xr.open_dataset(u_clim_file)
vclim_ds = xr.open_dataset(v_clim_file)

qclim = qclim_ds["var133"]
uclim = uclim_ds["var131"]
vclim = vclim_ds["var132"]

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

    # Match daily climatology by day of year
    doy = q["time"].dt.dayofyear

    q_clim_daily = qclim.sel(dayofyear=doy)
    u_clim_daily = uclim.sel(dayofyear=doy)
    v_clim_daily = vclim.sel(dayofyear=doy)

    # Replace dayofyear dimension with time
    q_clim_daily = q_clim_daily.assign_coords(time=q["time"]).swap_dims({"dayofyear": "time"})
    u_clim_daily = u_clim_daily.assign_coords(time=u["time"]).swap_dims({"dayofyear": "time"})
    v_clim_daily = v_clim_daily.assign_coords(time=v["time"]).swap_dims({"dayofyear": "time"})

    q_clim_daily = q_clim_daily.drop_vars("dayofyear", errors="ignore")
    u_clim_daily = u_clim_daily.drop_vars("dayofyear", errors="ignore")
    v_clim_daily = v_clim_daily.drop_vars("dayofyear", errors="ignore")

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
    # Decomposition
    # =========================

    # Dynamic term: q_clim * V'
    IVTx_dyn = vertical_integral(q_clim_daily * u_prime, dp, plev_name)
    IVTy_dyn = vertical_integral(q_clim_daily * v_prime, dp, plev_name)
    IVT_dyn = magnitude(IVTx_dyn, IVTy_dyn)

    # Thermodynamic/moisture term: q' * V_clim
    IVTx_thermo = vertical_integral(q_prime * u_clim_daily, dp, plev_name)
    IVTy_thermo = vertical_integral(q_prime * v_clim_daily, dp, plev_name)
    IVT_thermo = magnitude(IVTx_thermo, IVTy_thermo)

    # Nonlinear term: q' * V'
    IVTx_nonlinear = vertical_integral(q_prime * u_prime, dp, plev_name)
    IVTy_nonlinear = vertical_integral(q_prime * v_prime, dp, plev_name)
    IVT_nonlinear = magnitude(IVTx_nonlinear, IVTy_nonlinear)

    # =========================
    # Save only requested variables
    # =========================
    ds_out = xr.Dataset(
        {
            "IVT": IVT,
            "IVT_dyn": IVT_dyn,
            "IVT_thermo": IVT_thermo,
            "IVT_nonlinear": IVT_nonlinear,
        }
    )

    ds_out["IVT"].attrs["long_name"] = "Total integrated vapor transport"
    ds_out["IVT_dyn"].attrs["long_name"] = "Dynamic contribution to IVT anomaly: climatological moisture times anomalous wind"
    ds_out["IVT_thermo"].attrs["long_name"] = "Thermodynamic contribution to IVT anomaly: anomalous moisture times climatological wind"
    ds_out["IVT_nonlinear"].attrs["long_name"] = "Nonlinear contribution to IVT anomaly: anomalous moisture times anomalous wind"

    for var in ds_out.data_vars:
        ds_out[var].attrs["units"] = "kg m-1 s-1"

    out_file = f"{out_dir}/IVT_decomposition_{year}.nc"
    ds_out.to_netcdf(out_file)

    print(f"Saved: {out_file}")

    q_ds.close()
    u_ds.close()
    v_ds.close()
