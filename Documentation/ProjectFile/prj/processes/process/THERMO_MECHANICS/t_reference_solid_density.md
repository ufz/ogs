Reference solid density \f$\rho_{\rm SR0}\f$.

The real solid's density (the density of the solid grains) is expressed as:
\f[\rho_{\rm SR} = \frac{\rho_{\rm SR0}}{1 + 3 * \alpha_{\rm T} * \Delta T }\f]
under the assumption of mechanically incompressible solid.

The original Ansatz made is \f[\rho_{\rm SR} = \rho_{\rm SR0} {\rm exp} (-3
\alpha_{\rm T} \Delta T)\f] for large thermal strains. Once we linearize around
\f$\Delta T = 0\f$ and truncate the Taylor series after the first member, the
end result is the above formula \f$\rho_{\rm SR} = \rho_{\rm SR0} / ( 1 + 3
\alpha_{\rm T} \Delta T)\f$.  One can also arrive at the result by assuming a
linear volume strain from the start: \f$e = 3 \alpha_{\rm T} \Delta T\f$. The
key assumption is that the solid density is independent of pressure, _i.e._:
\f$\rho_{\rm SR} = \rho_{\rm SR}(T)\f$.
