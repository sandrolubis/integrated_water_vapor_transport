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
IVT_x = \frac{1}{g} \int_{p_{top}}^{p_{sfc}} q u \; dp
$$

and

$$
IVT_y = \frac{1}{g} \int_{p_{top}}^{p_{sfc}} q v \; dp
$$

Equivalently,

$$
IVT =
\left[
\left(
\frac{1}{g} \int_{p_{top}}^{p_{sfc}} q u \; dp
\right)^2
+
\left(
\frac{1}{g} \int_{p_{top}}^{p_{sfc}} q v \; dp
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
kg \; m^{-1} \; s^{-1}
$$

In discrete form, the IVT components can be approximated as:

$$
IVT_x = \frac{1}{g} \sum_{k} q_k u_k \Delta p_k
$$

$$
IVT_y = \frac{1}{g} \sum_{k} q_k v_k \Delta p_k
$$

where \(k\) represents each pressure level and \(\Delta p_k\) is the pressure-layer thickness.
