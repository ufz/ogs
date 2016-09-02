Sets up a convergence criterion that checks only the norm of the solution update
from one iteration to another.

Relative and/or absolute tolerance values can be specified, which are applied to
the norm of each component (i.e., each process variable, and, if there is a
vectorial process variable, to each vector component) of the solution vector
individually.

For each component at least one of the relative or absolute tolerance has to be
met in order to satisfy this convergence criterion, i.e., the following has to
hold (with errors \f$e\f$ and tolerances \f$t\f$):
\f[ e_{\mathrm{abs}} \le t_{\mathrm{abs}} \vee e_{\mathrm{rel}} \le t_{\mathrm{rel}} \f]
