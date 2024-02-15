This optional parameter indicates whether the porosity and transport porosity
should be initialized from the medium properties `porosity` and
`transport_porosity`, respectively.

This initialization is desired in almost all cases with one notable exception:
simulation restart when a non-constant porosity model is used,
e.g. with \ref ogs_file_param__prj__media__medium__properties__property__PorosityFromMassBalance "PorosityFromMassBalance"
or \ref ogs_file_param__prj__media__medium__properties__property__TransportPorosityFromMassBalance "TransportPorosityFromMassBalance".

In such a restart case simply set this parameter to `false`.
