## Integrated Water Vapor Transport (IVT) ##
by Sandro W. Lubis (PNNL)

This script computes integrated water vapor transport by integrating specific humidity (q) and horizontal wind components (u, v) over a specified vertical pressure range. The output can be used to study moisture flux patterns and their role in large-scale weather systems.

## How IVT is Calculated

The **Integrated Vapor Transport (IVT)** is computed as:

$\textbf{IVT} = \left[ \left( \frac{1}{g} \int_{p_{sfc}}^{p_{top}} qu , dp \right)^2 + \left( \frac{1}{g} \int_{p_{sfc}}^{p_{top}} qv , dp \right)^2 \right]^{\frac{1}{2}}$

where:

- **Zonal Component (IVTx):**
- 
  IVTx  =  \int_{p_{\text{top}}}^{p_{\text{sfc}}} q \cdot u \frac{dp}{g} 

- **Meridional Component (IVTy):**
- 
  IVTy =  \int_{p_{\text{top}}}^{p_{\text{sfc}}} q \cdot v \frac{dp}{g} 

Here:
- \( q \) = Specific humidity (kg/kg)
- \( u, v \) = Zonal and meridional wind components (m/s)
- \( dp \) = Pressure thickness of each level (must be in Pa)
- \( g \) = Gravitational acceleration (9.81 m/s²)
- \( p_{\text{sfc}} \) = Surface pressure
- \( p_{\text{top}} \) = Upper pressure limit


