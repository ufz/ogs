Under-relaxation (damping) coefficient for the Picard iteration.

Here \f$ x_k \f$ is the current iterate and \f$ g \f$ is the Picard fixed-point
map, i.e. \f$ g(x_k) \f$ is the solution of the linearized system
\f$ A(x_k)\,x = b(x_k) \f$ (the result of one linear solve); a converged iteration
is a fixed point \f$ x = g(x) \f$.

The default value 1.0 gives a non-damped Picard iteration
\f$ x_{k+1} = g(x_k) \f$. Values of the damping factor \f$ \beta \f$ in the range
(0, 1] provide under-relaxation for stabilization: the update is computed as
\f$ x_{k+1} = (1-\beta)\,x_k + \beta\,g(x_k) \f$.

Note on convergence checks: damping scales the increment by \f$ \beta \f$, i.e.
\f$ x_{k+1} - x_k = \beta\,(g(x_k) - x_k) \f$. A delta-x convergence criterion
therefore sees an increment reduced by \f$ \beta \f$ and may report convergence
prematurely for small \f$ \beta \f$. Tighten the delta-x tolerance accordingly,
or prefer a residual-based criterion, which is evaluated on the iterate before
damping is applied and is unaffected by the damping factor.
