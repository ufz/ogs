Template for the prefix which will be prepended to the output files. Allowed
template expressions are: `{:meshname}`, `{:timestep}`, `{:time}`,
`{:iteration}`.

Default value is `{:meshname}`.

Remark: The name of the output files will be constructed out of the prefix and
the suffix.

Furthermore, it is possible to specify the format of the expressions above. For
instance {:0>3timestep} results in a 3-digit output, if necessary with preceding
zeros. Also the time output can be formatted using the typical floating
point precision (for instance 0.4) and type (e, E, f, F, g, G) modifiers.
