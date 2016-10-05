Sets up a convergence criterion that checks only the norm of the residuum
vector.

A relative and/or absolute tolerance value can be specified, which is applied to
the norm of the whole residuum vector.

At least one of the relative or absolute tolerance has to be met in order to
satisfy this convergence criterion, i.e., the following has to hold (with errors
\f$e\f$ and tolerances \f$t\f$):
\f[ e_{\mathrm{abs}} \le t_{\mathrm{abs}} \vee e_{\mathrm{rel}} \le t_{\mathrm{rel}} \f]
