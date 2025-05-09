Template for the suffix which will be appended to the output files. Allowed
template expressions are: `{:meshname}`, `{:timestep}`, `{:time}`,
`{:iteration}`, `{:converged}`.

Default value is `_ts_{:timestep}_t_{:time}`.

Remark: The name of the output files will be constructed out of the prefix and
the suffix.

Furthermore, it is possible to specify the format of the expressions above. For
instance {:0>3timestep} results in a 3-digit output, if necessary with preceding
zeros. Also the time output can be formatted using the typical floating
point precision (for instance 0.4) and type (e, E, f, F, g, G) modifiers.

If the non-linear solver norms are not met the `{:converged}` expression is replaced with
`_not_converged` and removed otherwise.
