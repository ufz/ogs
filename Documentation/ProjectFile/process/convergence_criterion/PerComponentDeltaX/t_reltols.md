Relative tolerance values being applied.

In this list of values there must be one value for each process variable (and if
there are vectorial process variables, for each vector component of that process
variable). The order of values must match the order of process variables as
defined by the specific process.

If for a certain component no relative tolerance shall be applied, zero can be
given as the respective tolerance value.

For the computation of the relative error for each component, see \ref
ogs_file_param__process__convergence_criterion__DeltaX__reltol.
