## Integrated Water Vapor Transport (IVT) ##
by Sandro W. Lubis (PNNL)

This script computes integrated water vapor transport by integrating specific humidity (q) and horizontal wind components (u, v) over a specified vertical pressure range. The output can be used to study moisture flux patterns and their role in large-scale weather systems.

## How IVT is Calculated

The **Integrated Vapor Transport (IVT)** is computed as:

<p align="center"> 
    \( IVT = \sqrt{IVTx^2 + IVTy^2} \)
</p>

where:

- **Zonal Component (IVTx)**:
  <p align="center"> 
      \( IVTx = \int_{p_{\text{top}}}^{p_{\text{sfc}}} q \cdot u \frac{dp}{g} \)
  </p>

- **Meridional Component (IVTy)**:
  <p align="center"> 
      \( IVTy = \int_{p_{\text{top}}}^{p_{\text{sfc}}} q \cdot v \frac{dp}{g} \)
  </p>

Here:
- \( q \) = Specific humidity (kg/kg)
- \( u, v \) = Zonal and meridional wind components (m/s)
- \( dp \) = Pressure thickness of each level
- \( g \) = Gravitational acceleration (9.81 m/s²)
- \( p_{\text{sfc}} \) = Surface pressure
- \( p_{\text{top}} \) = Upper pressure limit (300 hPa in this case)


## IVT Calculation

The **Integrated Vapor Transport (IVT)** is computed as:

<div align="center">
  \( IVT = \sqrt{IVTx^2 + IVTy^2} \)
</div>

where:

- **Zonal Component (IVTx):**
  <div align="center"> 
      \( IVTx = \int_{p_{\text{top}}}^{p_{\text{sfc}}} q \cdot u \frac{dp}{g} \)
  </div>

- **Meridional Component (IVTy):**
  <div align="center"> 
      \( IVTy = \int_{p_{\text{top}}}^{p_{\text{sfc}}} q \cdot v \frac{dp}{g} \)
  </div>

### Enabling MathJax on GitHub Pages
To enable MathJax on GitHub Pages, add this to your `_config.yml` file:

```yml
markdown: kramdown
kramdown:
  math_engine: mathjax

