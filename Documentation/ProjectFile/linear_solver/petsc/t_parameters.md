Sets PETSc options from the given parameter string.

The default value is <tt>"-ksp_type cg -pc_type bjacobi -ksp_rtol 1e-16 -ksp_max_it 10000"</tt>

\attention
PETSc can also receive options from the commandline arguments. However, the
values from the prj file are applied directly to the linear solver. If
specifying both, make sure that you check which options take precedence.
