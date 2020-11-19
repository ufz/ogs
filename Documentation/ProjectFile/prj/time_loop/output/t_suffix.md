Template for the suffix which will be appended to the output files. Allowed
template expressions are:
- {:meshname}
- {:timestep}
- {:time}.

Default value: _ts_{:timestep}_t_{:time}

Remark: The name of the result files will be constructed out of the prefix and
the suffix.

Furthermore, it is possible to specifiy the format of the expressions above. For
instance {:0>3timestep} results in a 3-digit output, if necessary with preceding
zeros. Also the time output can be formated using the typical floating
point precision (for instance 0.4) and type (e, E, f, F, g, G) modifiers.
