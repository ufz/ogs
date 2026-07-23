Depth of Anderson acceleration for the Picard iteration, i.e. the number of
previous iterates \f$ m \f$ retained to accelerate convergence.

Here \f$ x_k \f$ is the current iterate and \f$ g \f$ is the Picard fixed-point
map, i.e. \f$ g(x_k) \f$ is the solution of the linearized system
\f$ A(x_k)\,x = b(x_k) \f$ (the result of one linear solve); a converged
iteration is a fixed-point \f$ x = g(x) \f$.

For \f$ m \ge 2 \f$ the next iterate is a mixture of the last \f$ m \f$ (damped)
steps whose weights minimise the residual \f$ \|\sum_i \theta_i f_i\| \f$
subject to \f$ \sum_i \theta_i = 1 \f$, which typically reduces the number of
iterations for slowly converging fixed-point problems.

The values \f$ m = 0 \f$ and \f$ m = 1 \f$ are accepted but have no effect:
mixing a single stored step yields \f$ \theta = (1) \f$ and hence reproduces the
plain Picard update \f$ x_{k+1} = g(x_k) \f$ (subject to the `damping` factor).
A warning is emitted in that case; omit the `anderson` subtree to disable the
acceleration explicitly.

Whenever the stored steps become (nearly) linearly dependent the mixing weights
are not trustworthy. Such an iteration falls back to the plain Picard step and
says so with an `info` message, so the acceleration degrades gracefully instead
of amplifying rounding error.

The value must be non-negative. Anderson acceleration is not compatible with
linear equation systems (a single Picard step already yields the exact
solution); omit the `anderson` subtree for those.
