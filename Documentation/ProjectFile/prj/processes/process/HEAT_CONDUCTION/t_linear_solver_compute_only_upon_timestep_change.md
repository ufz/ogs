This flag enables an optimization for the `HEAT_CONDUCTION` process:

the linear solver will only do the `compute()`<sup>\*</sup> step if the time step size
changes.

This flag is a further optimization on top of the
[\<linear\>](@ref ogs_file_param__prj__processes__process__HEAT_CONDUCTION__linear)
flag.
Thus, the requirements of `<linear>` apply to this flag, too!

\attention
This is an expert option. It comes with a number of further **requirements above
those of `<linear>`**. These are:

- The linear solver used to solve the process equations must be exclusively used
  for a single process (or for a single `process_id` in the case of staggered
  coupling). Otherwise the equation system the linear solver has to solve would
  not only change upon time step size change, but also when switching processes.
- Both, the process equations and the contributions from source terms and
  boundary conditions must be **solution and time independent**!
  In particular, the matrices presented to the linear solver must only change
  upon time step size change.
  Note that some boundary conditions (Robin, Python, ...) might contribute not only
  to the right-hand side of the equation system, but also to the matrix.
  The contributions to the matrix must be solution and time independent!
  Also, switching Dirichlet BCs on and off for certain d.o.f.s (e.g., via Python
  boundary conditions) is not possible!

\attention
OGS does not detect whether any or all of the above conditions are met.
It is the sole responsibility of the user to ensure the correct usage of this
expert option!
The user is **strongly advised to check the correct use of this option** by
comparing results to simulation runs that do not use this option!

<sup>\*</sup>: `compute()` computes the LU decomposition for a direct linear
solver, or computes the preconditioner for an iterative linear solver.
