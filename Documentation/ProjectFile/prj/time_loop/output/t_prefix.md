Template for the prefix which will be prepended to the output files. Allowed
template expressions are:
- {:meshname}
- {:timestep}
- {:time}.

Default value: {:meshname}

Remark: The name of the pvd-file will be constructed out of the prefix.

Furthermore, it is possible to specifiy the format of the expressions above. For
instance {:0>3timestep} results in a 3-digit output, if necessary with preceding
zeros. Also the time output can be formated using the typical floating
point precision (for instance 0.4) and type (e, E, f, F, g, G) modifiers.
