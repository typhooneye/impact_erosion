# Atmospheric Erosion

## 1. Planetesimal Impacts

### Assumptions

1. The impactor has constant density.
2. The impactor is spherical.
3. For the target planet, the radius-mass scaling relation is assumed to be
   \[
   \left(\frac{R}{R_E}\right) = \left(\frac{M}{M_E}\right)^{0.27}
   \]
   (reference needed).
4. The atmospheric scale height is calculated from the surface gravity.

### 1.1 Planetesimal Impact Properties

| Property | Symbol |
| :--- | :---: |
| Size | \(D\) |
| Density | \(\rho_{\mathrm{imp}}\) |
| Velocity | \(v_{\mathrm{imp}}\) |
| Mass | \(m_{\mathrm{imp}} = \frac{\pi}{6}\rho_{\mathrm{imp}}D^3\) |

### 1.2 Target Properties

| Property | Symbol |
| :--- | :---: |
| Radius | \(R_p\) |
| Mass | \(M_p\) |
| Density | \(\rho_p\) |
| Surface gravity | \(g\) |
| Mean molecular weight | \(\mu\) |
| Atmospheric density | \(\rho_{\mathrm{atm}}\) |
| Atmospheric temperature | \(T_{\mathrm{atm}}\) |
| Escape velocity | \(v_{\mathrm{esc}} = \sqrt{\frac{2GM}{R}}\) |
| Atmospheric scale height | \(H = \frac{k_B T}{\mu m_H g}\) |
| Atmospheric mass | \(m_{\mathrm{atm}} = 4\pi R_p^2 H \rho_{\mathrm{atm}}\) |

**Open question:** When does atmospheric escape begin? Can this atmospheric mass estimate be used consistently for that stage?

### 1.3 Size Distribution

\[
\frac{dN}{dD} = kD^{-2}
\]

To solve for \(k\):

\[
\int_{D_{\min}}^{D_{\max}} \frac{dN}{dD}\frac{\pi}{6}\rho D^3 \, dD = M_{\mathrm{tot}}
\]

\[
\int_{D_{\min}}^{D_{\max}} \frac{dN}{dD} D^3 \, dD = \frac{6M_{\mathrm{tot}}}{\pi \rho}
\]

\[
\int_{D_{\min}}^{D_{\max}} kD \, dD = \frac{6M_{\mathrm{tot}}}{\pi \rho}
\]

\[
k = \frac{12M_{\mathrm{tot}}}{\pi \rho \left(D_{\max}^2 - D_{\min}^2\right)}
\]

Then,

\[
N = \int_{D_{\min}}^{D_{\max}} kD^{-2} \, dD = k\left(\frac{1}{D_{\min}} - \frac{1}{D_{\max}}\right)
\]

## 2. Giant Impacts

### Assumptions

1. A scaled radius-mass relation is assumed for giant impactors:
   \[
   \left(\frac{R}{R_E}\right) = \left(\frac{M}{M_E}\right)^{0.27}
   \]
   (reference needed).

### 2.1 Giant Impact Properties

| Property | Symbol |
| :--- | :---: |
| Radius | \(R_i\) |
| Bulk density | \(\rho_{\mathrm{imp}}\) |
| Velocity | \(v_{\mathrm{imp}}\) |
| Mass | \(m_{\mathrm{imp}} = \frac{\pi}{6}\rho_{\mathrm{imp}}D^3\) |
| Impact parameter | \(b = \sin(\beta)\) |
