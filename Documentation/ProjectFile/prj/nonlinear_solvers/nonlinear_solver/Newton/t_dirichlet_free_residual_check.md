If set to `true`, the residual passed to the convergence criterion excludes the
Dirichlet degrees of freedom. Default: `false`.

The checked vector is \f$ P\,(r - J V) \f$, where \f$ V_k = x_k - g_k \f$ holds
the deviations of the current solution from the Dirichlet values \f$ g_k \f$ at
the Dirichlet degrees of freedom \f$ k \f$ and is zero elsewhere, and \f$ P \f$ sets the
entries at the Dirichlet degrees of freedom to zero. This is the residual linearized at the
state fulfilling the Dirichlet BCs, restricted to the free degrees of freedom.

Compared with the default, the norm contains neither the reaction forces nor
the Dirichlet rows of the linear system, whose scaling is arbitrary and
differs between the linear algebra backends. In an iteration in which the
Dirichlet values changed, e.g., the first one of a time step, the checked
residual still contains the out-of-balance forces caused by that change.
Therefore, it cannot be satisfied before the Dirichlet BCs are imposed.

Relative residual tolerances refer to a different \f$ |r_0| \f$ than with the
default setting.
