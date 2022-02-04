\copydoc NumLib::NonlinearSolver<NonlinearSolverTag::Newton>::_compensate_non_equilibrium_initial_residuum

If this option tag is not given or set to false for all variables, the calculation and the use of
 the initial non equilibrium initial residuum is skipped. Otherwise, the calculation of the initial
 non equilibrium initial residuum is triggered, and the portions of the calculated initial residual
 vector, which are associated with the variables without this option enabled are set to zero.

With this equation variable option, one can flexibly set which portions of the equations can
 have initial non_equilibrium compensation. For example, in a numerical simulation of
 coupled THM processes in porous medium, if you want to avoid the initial stress
  influence the displacement result, you can simply set this option to true only
 for the variable of displacement in the variable configuration in the project file.
