The relative tolerance being applied to the norm of the whole residuum vector.

The relative error \f$e\f$ is computed as
\f[e = \frac{||r||}{||r_0||+\epsilon},\f]
where \f$r\f$ and \f$r_0\f$ are the residuum and the residuum from computed in
the first iteration, respectively, and \f$\epsilon\f$ is the machine epsilon
guaranteeing that the denominator is always nonzero.
