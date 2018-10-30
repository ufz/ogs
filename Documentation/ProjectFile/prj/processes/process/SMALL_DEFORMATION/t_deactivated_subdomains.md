Set of material IDs. The subdomains associated to those material IDs are
 excluded from the local assembly.

It this option is used, please make sure that there is no boundary condition
 at any node inside these subdomains, and only the iterative linear solver is
 utilized. If there are boundary conditions  at the nodes inside these
 subdomains, the program will terminate with linear solver failure.
