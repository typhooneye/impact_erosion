# Method: Atmospheric Erosion

During planet formation, atmospheric mass loss can occur due to both giant impacts from embryos and planetesimal accretion. We adopt the prescriptions for planetesimal-sized impacts (1 m to 100 km; \(<10^{18}\) kg) from Shuvalov (2009) and Schlichting et al. (2015), and for embryo-sized impacts (0.005 to 0.48 Earth masses) from Kegerreis (2020).

## Planetesimal-Sized Impacts

For planetesimal-sized impacts, the N-body simulations account only for the cumulative effects of a series of planetesimals with total mass \(M_{\mathrm{tot}}\). To resolve the impacts, the impactors (\(N\)) are assumed to follow a power-law size distribution with a differential power-law index of \(q = 2\), as inferred from the lunar cratering record by Hartmann (1964), or a differential power-law index of \(q = 3.5\), which is the typical value for a self-similar cascade of fragments (Dohnanyi 1969; Neukum et al. 2001):

\[
\frac{dN}{dD} \propto kD^{-2},
\]

where \(N\) is the number of impactors and \(D\) is the diameter of the impactor. The largest impactor has \(D_{\max} = 300\,\mathrm{km}\), and the smallest has \(D_{\min} = 1\,\mathrm{m}\) (Schlichting and Mukhopadhyay 2018). We assume a planetesimal density of

\[
\rho_i = 3\,\mathrm{g\,cm^{-3}},
\]

and the mass of a planetesimal of size \(D_i\) is

\[
m_i = \frac{\pi \rho_i D_i^3}{6}.
\]

We convert the size distribution to a mass distribution and integrate it to match the total mass. The total number of impactors for a given total mass \(M_{\mathrm{tot}}\) is then

\[
N = \frac{12M_{\mathrm{tot}}}{\pi \rho_i \left(D_{\max}D_{\min}(D_{\max}+D_{\min})\right)}.
\]

We then generate a random sample following the given size distribution with sample size \(N\).

During planetesimal accretion, atmospheric mass loss can be induced by crater formation at the impact site if the impactor reaches the surface intact, or by aerial burst and/or impactor fragmentation during its passage through the atmosphere. Assuming an isothermal atmosphere, Shuvalov et al. (2014) proposed two key parameters to distinguish these effects: \(H_{1/2}\), the altitude at which the projectile velocity is reduced by half and where aerial bursts are assumed to occur, and \(H_{\mathrm{fr}}\), the fragmentation altitude for disrupting projectiles:

\[
H_{1/2} = -A_{1/2}H \ln\!\left(\frac{B_{1/2}}{2H}D\sqrt{\frac{\rho_i}{\rho_0}}\right),
\]

\[
H_{\mathrm{fr}} = -A_{\mathrm{fr}}H \ln\!\left(\frac{B_{\mathrm{fr}}}{2H}D\sqrt{\frac{\rho_i}{\rho_0}}\right).
\]

Here, \(H\) is the atmospheric scale height, \(D\) is the impactor size, \(\rho_i\) is the impactor density, and \(\rho_0\) is the atmospheric density at the base of the atmosphere. The scale factors are obtained by fitting numerical simulations. Here we adopt either the Shuvalov et al. (2014) values

\[
A_{1/2} = A_{\mathrm{fr}} = 1.4,\quad B_{1/2} = 5.5,\quad B_{\mathrm{fr}} = 3.5,
\]

or the modified values from Sinclair and Wyatt (2021),

\[
A_{1/2} = A_{\mathrm{fr}} = 1.4,\quad B_{1/2} = 1.24,\quad B_{\mathrm{fr}} = 1.97.
\]

The impact outcome is classified as follows:

- If \(H_{1/2} < 0\) and \(H_{\mathrm{fr}} < 0\), the impactor reaches the surface intact.
- If \(H_{1/2} < 0\) and \(H_{\mathrm{fr}} > 0\), the impactor fragments in the atmosphere.
- If \(H_{1/2} > 0\) and \(H_{\mathrm{fr}} > 0\), the impactor undergoes an aerial burst.

### Cratering Impacts

For cratering impacts (\(H_{1/2}, H_{\mathrm{fr}} < 0\)), atmospheric loss is parameterized using the prescription from Shuvalov (2009), based on angle-averaged numerical simulations. The prescription is written in terms of the dimensionless erosional efficiency \(\xi\), which depends on the impactor properties (size \(D\), density \(\rho_i\), and velocity \(v_{\mathrm{imp}}\)) and target properties (escape velocity \(v_{\mathrm{esc}}\), surface density \(\rho_p\), atmospheric scale height \(H\), and atmospheric density at the base of the atmosphere \(\rho_0\)):

\[
\xi = \left(\frac{D}{H}\right)^3
\left[\left(\frac{v_{\mathrm{imp}}}{v_{\mathrm{esc}}}\right)^2 - 1\right]
\left[\frac{\rho_i \rho_p}{\rho_0(\rho_i + \rho_p)}\right].
\]

The scale height of the atmosphere is

\[
H = \frac{k_B T}{\mu g},
\]

where \(k_B\) is the Boltzmann constant, \(T\) is the surface temperature, \(\mu\) is the mean molecular weight, and \(g\) is the gravitational acceleration as a function of planetary radius and mass. We assume an isothermal atmosphere, and the adopted values are shown in Table X. This isothermal assumption is necessary in order to apply the Shuvalov (2009) prescription.

The atmospheric mass loss due to a single impact is

\[
\left(\frac{m_{\mathrm{atm\,loss}}}{m_i}\right)_{\mathrm{cr}}
=
\left[\left(\frac{v_{\mathrm{imp}}}{v_{\mathrm{esc}}}\right)^2 - 1\right]\chi_a(\xi),
\]

where

\[
\log_{10}\!\left(\chi_a(\xi)\right)
=
-6.375 + 5.239\log_{10}\xi - 2.121(\log_{10}\xi)^2 + 0.397(\log_{10}\xi)^3 - 0.037(\log_{10}\xi)^4 + 0.0013(\log_{10}\xi)^5.
\]

Following Sinclair et al. (2020) and Sinclair and Wyatt (2021), we modify this expression for \(\xi > 10^6\):

\[
\log_{10}\!\left(\chi_a(\xi > 10^6)\right) = -0.6438\xi + 0.4746.
\]

We also adopt the polar-cap limit from Schlichting et al. (2015), such that the atmospheric loss from a single impact cannot exceed the mass of the atmosphere contained in the polar cap:

\[
m_{\max} = 2\pi \rho_0 H^2 R_p,
\]

where \(R_p\) is the radius of the planet.

### Aerial Burst

For aerial bursts (\(H_{1/2} > 0\), \(H_{\mathrm{fr}} > 0\)), the atmospheric loss is given by Shuvalov et al. (2014):

\[
\left(\frac{m_{\mathrm{atm\,loss}}}{m_i}\right)_{\mathrm{ab}}
=
5.5\times10^{-3}(v_{\mathrm{imp}} - 15 - 0.2H)^2\sqrt{\rho_i}\left(\frac{D}{H}\right).
\]

Here, \(v_*\) is the minimum entry velocity required to produce atmospheric erosion. Note that \(\rho_i\) is in units of \(\mathrm{g\,cm^{-3}}\), velocities are in \(\mathrm{km\,s^{-1}}\), and distance is in km. The threshold velocity scales with atmospheric scale height as

\[
v_* = 15 + 0.2H.
\]

### Fragmentation

For fragmentation (\(H_{1/2} < 0\), \(H_{\mathrm{fr}} > 0\)), the mass loss is interpolated between the cratering and aerial-burst regimes:

\[
\ln\!\left(\frac{m_{\mathrm{atm\,loss}}}{m_i}\right)_{\mathrm{fr}}
=
\ln\!\left(\frac{m_{\mathrm{atm\,loss}}}{m_i}\right)_{\mathrm{cr}}
+
\left[
\ln\!\left(\frac{m_{\mathrm{atm\,loss}}}{m_i}\right)_{\mathrm{ab}}
-
\ln\!\left(\frac{m_{\mathrm{atm\,loss}}}{m_i}\right)_{\mathrm{cr}}
\right]
\frac{H_{\mathrm{fr}}}{H_{\mathrm{fr}} - H_{1/2}}.
\]

Note from the original draft: in Sinclair's thesis, there is a more complicated modification to this expression. It is still unclear whether that modification is needed.

## Giant Impacts

The above prescription applies to impacts during planetesimal accretion, while giant impacts can erode the atmosphere globally by driving a shock wave through the planet. We follow the scaling law obtained from numerical simulations that consider different impactor-target mass ratios, densities, and impact angles (Kegerreis et al. 2020).

The fraction of eroded atmosphere, \(X\), is

\[
\frac{m_{\mathrm{atm\,loss}}}{m_i}\frac{M_a + M_p}{m_i}
=
0.64
\left[
\left(\frac{v_c}{v_{\mathrm{esc}}}\right)^2
\left(\frac{m_i}{M_a + M_p}\right)^{1/2}
\left(\frac{\rho_i}{\rho_p}\right)^{1/2}
f_M(b)
\right]^{0.65}.
\]

Here, \(M_a\) is the atmospheric mass, \(M_p\) is the planetary mass, and \(f_M(b)\) is the fractional interacting mass defined as

\[
f_M(b) \equiv \frac{\rho_p V_p^{\mathrm{cap}} + \rho_i V_i^{\mathrm{cap}}}{\rho_p V_p + \rho_i V_i}.
\]

Here, \(V_p\) and \(V_i\) are the total volumes of the planet and impactor, and \(V_p^{\mathrm{cap}}\) and \(V_i^{\mathrm{cap}}\) are the cap volumes defined by the geometry at contact. Specifically, \(V_p^{\mathrm{cap}}\) refers to the volume of the planet above the lowest point of the impactor at contact, and \(V_i^{\mathrm{cap}}\) refers to the volume of the impactor below the highest point of the planet.

The height of both caps is

\[
d_{\mathrm{cap}} = \left(R_p + \frac{D}{2}\right)(1-b),
\]

which gives

\[
V_{p/i}^{\mathrm{cap}} = \frac{\pi d_{\mathrm{cap}}^2}{3}(3R_{p/i} - d_{\mathrm{cap}}).
\]

The heat released by giant impacts also causes thermal expansion of the atmosphere, which may drive atmospheric escape through a Parker-type wind, as examined by Biersteker and Schlichting (2021). They show that only if the atmosphere's mean molecular weight is larger than 16 can an impactor drive thermal loss. The mean molecular weight in all of our simulations is above `XXX`, so thermal loss should be negligible.

## XUV-Driven Loss

TBD.

## References

- Biersteker J. B., Schlichting H. E., 2021, MNRAS, 501, 587
- Dohnanyi, J. S. 1969, J. Geophys. Res., 74, 2531
- Kegerreis, J. A., Eke, V. R., Catling, D. C., et al. 2020a, ApJ, 901, L31
- Hartmann, W. K. 1964, "On the Distribution of Lunar Crater Diameters", Comm. Lunar Planet. Lab. 2, 197-203
- Neukum, G., Ivanov, B. A., Hartmann, W. K. 2001, "Cratering records in the inner solar system in relation to the lunar reference system", Space Sci. Rev. 96, 55-86
- Shuvalov, V., Kührt, E., de Niem, D., Wünnemann, K. 2014, Planet. Space Sci., 98, 120
- Schlichting, H. E., Mukhopadhyay, S. 2018, "Atmosphere Impact Losses", Space Sci. Rev. 214, 34. https://doi.org/10.1007/s11214-018-0471-z
- Sinclair, C. and Wyatt, M. 2021, EGU General Assembly Conference Abstracts. doi:10.5194/egusphere-egu21-2628
- Sinclair, C. A., Wyatt, M. C., Morbidelli, A., Nesvorný, D. 2020, MNRAS, 499, 5334
- Sinclair, C. A., Wyatt, M. C. 2021, MNRAS, 509, 345
