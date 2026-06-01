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

The IVT can also be decomposed into contributions from changes in moisture, changes in winds, and their nonlinear interaction.

Let specific humidity and horizontal winds be separated into climatological and anomalous components:

$$
q = \overline{q} + q'
$$

$$
u = \overline{u} + u'
$$

$$
v = \overline{v} + v'
$$

where the overbar denotes the daily climatology and the prime denotes the anomaly relative to that climatology.

The zonal IVT component is

$$
IVT_x = \frac{1}{g} \int_{p_{top}}^{p_{sfc}} q u \; dp
$$

Substituting \(q = \overline{q} + q'\) and \(u = \overline{u} + u'\):

$$
IVT_x =
\frac{1}{g} \int_{p_{top}}^{p_{sfc}}
(\overline{q} + q')(\overline{u} + u') \; dp
$$

Expanding the terms gives

$$
IVT_x =
\frac{1}{g} \int_{p_{top}}^{p_{sfc}}
\overline{q}\overline{u} \; dp
+
\frac{1}{g} \int_{p_{top}}^{p_{sfc}}
\overline{q}u' \; dp
+
\frac{1}{g} \int_{p_{top}}^{p_{sfc}}
q'\overline{u} \; dp
+
\frac{1}{g} \int_{p_{top}}^{p_{sfc}}
q'u' \; dp
$$

Therefore, the IVT anomaly can be decomposed as

$$
IVT_x' =
IVT_{x,dyn}
+
IVT_{x,thermo}
+
IVT_{x,nonlinear}
$$

where

$$
IVT_{x,dyn}
=
\frac{1}{g} \int_{p_{top}}^{p_{sfc}}
\overline{q}u' \; dp
$$

$$
IVT_{x,thermo}
=
\frac{1}{g} \int_{p_{top}}^{p_{sfc}}
q'\overline{u} \; dp
$$

$$
IVT_{x,nonlinear}
=
\frac{1}{g} \int_{p_{top}}^{p_{sfc}}
q'u' \; dp
$$

Similarly, for the meridional component:

$$
IVT_y' =
IVT_{y,dyn}
+
IVT_{y,thermo}
+
IVT_{y,nonlinear}
$$

where

$$
IVT_{y,dyn}
=
\frac{1}{g} \int_{p_{top}}^{p_{sfc}}
\overline{q}v' \; dp
$$

$$
IVT_{y,thermo}
=
\frac{1}{g} \int_{p_{top}}^{p_{sfc}}
q'\overline{v} \; dp
$$

$$
IVT_{y,nonlinear}
=
\frac{1}{g} \int_{p_{top}}^{p_{sfc}}
q'v' \; dp
$$

The total anomaly contribution from each term is then calculated as the vector magnitude:

$$
IVT_{dyn}
=
\sqrt{
IVT_{x,dyn}^{2}
+
IVT_{y,dyn}^{2}
}
$$

$$
IVT_{thermo}
=
\sqrt{
IVT_{x,thermo}^{2}
+
IVT_{y,thermo}^{2}
}
$$

$$
IVT_{nonlinear}
=
\sqrt{
IVT_{x,nonlinear}^{2}
+
IVT_{y,nonlinear}^{2}
}
$$

where:

- **Dynamic term**: contribution from anomalous winds acting on climatological moisture  
- **Thermodynamic term**: contribution from anomalous moisture transported by climatological winds  
- **Nonlinear term**: contribution from the interaction between anomalous moisture and anomalous winds  

In the script, these terms are computed as:

```python
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
