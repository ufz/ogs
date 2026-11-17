Computes the normal flux across a surface mesh after each time step and
writes it as a cell property of that mesh.

For every surface element \f$\Gamma_e\f$ the stored value is the
element-averaged normal flux
\f[
  \bar{q}_e = \frac{1}{|\Gamma_e|} \int_{\Gamma_e} \mathbf{q} \cdot \mathbf{n}
  \,\mathrm{d}\Gamma,
\f]
with \f$\mathbf{n}\f$ the unit normal derived from the node order of the
surface element. For surface elements that keep the node order of the bulk
element faces, e.g. meshes extracted with `ExtractBoundary`, \f$\mathbf{n}\f$
points into the domain: inflow is positive, outflow negative. The flux
\f$\mathbf{q}\f$ is evaluated in the adjacent bulk element and its meaning and
unit depend on the process:

- LIQUID_FLOW: the Darcy flux including the gravity term in m/s for
  `equation_balance_type` `volume`, the mass flux \f$\rho\mathbf{q}\f$ in
  kg/(m^2 s) for `mass`;
- HT: the Darcy flux including the gravity term in m/s;
- ComponentTransport: the mass flux \f$\rho\mathbf{q}\f$ in kg/(m^2 s);
- SteadyStateDiffusion: the diffusive flux \f$-D\nabla u\f$.

Multiply by the element area to obtain the flow rate through the element.
