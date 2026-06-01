# Integrated Water Vapor Transport (IVT)

**Author:** Sandro W. Lubis, PNNL

This script calculates **Integrated Water Vapor Transport (IVT)** by vertically integrating specific humidity (`q`) and horizontal wind components (`u`, `v`) over a selected pressure layer. The resulting IVT fields can be used to analyze moisture transport, atmospheric rivers, and the role of moisture fluxes in large-scale weather systems.

## How IVT Is Calculated

The magnitude of IVT is computed from its zonal and meridional components:

$$
IVT = \sqrt{IVT_x^2 + IVT_y^2}
$$

where

$$
IVT_x = \frac{1}{g} \int_{p_{top}}^{p_{sfc}} q u \, dp
$$

and

$$
IVT_y = \frac{1}{g} \int_{p_{top}}^{p_{sfc}} q v \, dp
$$

Equivalently,

$$
\mathbf{IVT} =
\left[
\left(
\frac{1}{g} \int_{p_{top}}^{p_{sfc}} q u \, dp
\right)^2
+
\left(
\frac{1}{g} \int_{p_{top}}^{p_{sfc}} q v \, dp
\right)^2
\right]^{1/2}
$$

where:

- `q` is specific humidity in kg kg\(^{-1}\)
- `u` is the zonal wind component in m s\(^{-1}\)
- `v` is the meridional wind component in m s\(^{-1}\)
- `dp` is the pressure thickness of each layer in Pa
- `g` is gravitational acceleration, approximately 9.81 m s\(^{-2}\)
- \(p_{sfc}\) is the surface pressure
- \(p_{top}\) is the upper pressure limit of the integration layer

The units of IVT are:

$$
kg \, m^{-1} \, s^{-1}
$$



## IVT Anomaly Decomposition

This script decomposes IVT into dynamic, thermodynamic, and nonlinear contributions by separating moisture and wind into climatological and anomalous components.

Specific humidity and winds are written as:

$$q = \overline{q} + q'$$

$$u = \overline{u} + u'$$

$$v = \overline{v} + v'$$

where:

- $\overline{q}$ is the daily climatological specific humidity
- $q'$ is the specific humidity anomaly
- $\overline{u}$ and $\overline{v}$ are the daily climatological zonal and meridional winds
- $u'$ and $v'$ are the anomalous zonal and meridional winds

The zonal IVT component is:

$$IVT_x = \frac{1}{g} \int_{p_{top}}^{p_{sfc}} q u, dp$$

Substituting the climatological and anomalous components gives:

$$IVT_x = \frac{1}{g} \int_{p_{top}}^{p_{sfc}} (\overline{q} + q')(\overline{u} + u'), dp$$

After expansion:

$$IVT_x = \frac{1}{g} \int_{p_{top}}^{p_{sfc}} \overline{q}\overline{u}, dp + \frac{1}{g} \int_{p_{top}}^{p_{sfc}} \overline{q}u', dp + \frac{1}{g} \int_{p_{top}}^{p_{sfc}} q'\overline{u}, dp + \frac{1}{g} \int_{p_{top}}^{p_{sfc}} q'u', dp$$

Therefore, the zonal IVT anomaly can be decomposed as:

$$IVT_x' = IVT_{x,dyn} + IVT_{x,thermo} + IVT_{x,nonlinear}$$

where:

$$IVT_{x,dyn} = \frac{1}{g} \int_{p_{top}}^{p_{sfc}} \overline{q}u', dp$$

$$IVT_{x,thermo} = \frac{1}{g} \int_{p_{top}}^{p_{sfc}} q'\overline{u}, dp$$

$$IVT_{x,nonlinear} = \frac{1}{g} \int_{p_{top}}^{p_{sfc}} q'u', dp$$

Similarly, the meridional IVT anomaly is:

$$IVT_y' = IVT_{y,dyn} + IVT_{y,thermo} + IVT_{y,nonlinear}$$

where:

$$IVT_{y,dyn} = \frac{1}{g} \int_{p_{top}}^{p_{sfc}} \overline{q}v', dp$$

$$IVT_{y,thermo} = \frac{1}{g} \int_{p_{top}}^{p_{sfc}} q'\overline{v}, dp$$

$$IVT_{y,nonlinear} = \frac{1}{g} \int_{p_{top}}^{p_{sfc}} q'v', dp$$

The magnitude of each contribution is calculated as:

$$IVT_{dyn} = \sqrt{IVT_{x,dyn}^{2} + IVT_{y,dyn}^{2}}$$

$$IVT_{thermo} = \sqrt{IVT_{x,thermo}^{2} + IVT_{y,thermo}^{2}}$$

$$IVT_{nonlinear} = \sqrt{IVT_{x,nonlinear}^{2} + IVT_{y,nonlinear}^{2}}$$

where:

- **Dynamic contribution**: anomalous winds acting on climatological moisture
- **Thermodynamic contribution**: anomalous moisture transported by climatological winds
- **Nonlinear contribution**: interaction between anomalous moisture and anomalous winds

The output variables are:

- `IVT`: total integrated vapor transport
- `IVT_dyn`: dynamic contribution
- `IVT_thermo`: thermodynamic or moisture contribution
- `IVT_nonlinear`: nonlinear interaction contribution

All output variables have units of:

$$kg m^{-1} s^{-1}$$

This decomposition helps identify whether IVT changes are mainly caused by wind anomalies, moisture anomalies, or their combined effect.
The decomposition is exact for the vector components `IVT_x` and `IVT_y`. The scalar magnitudes of each term are useful for diagnosis, but they do not necessarily add exactly to the total IVT anomaly because IVT magnitude is calculated using a square root.
