The relative tolerance being applied to the norm of the whole solution vector.

The relative error \f$e\f$ is computed as
\f[e = \frac{||\Delta x||}{||x||+\epsilon},\f]
where \f$\Delta x\f$ and \f$x\f$ are the current solution update and the
new solution from the current iteration, respectively, and \f$\epsilon\f$ is the
machine epsilon guaranteeing that the denominator is always nonzero.
