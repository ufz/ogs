Set of material indices, the subdomains by which are excluded from local assembly.

It this option is used, please make sure that there is no boundary condition
 at any node inside these subdomains, and only the iterative linear solver is
 utilized. If there are boundary conditions  at the nodes inside these
 subdomains, the program will exit by the failure of the iterative linear.
