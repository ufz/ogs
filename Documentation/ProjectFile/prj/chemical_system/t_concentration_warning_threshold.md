Optional. A negative (or zero) reporting threshold for clamped concentrations.
All negative input concentrations are clamped to zero, because PHREEQC rejects
negatives; this parameter only sets the value below which the clamping is
reported, it does not gate the clamping itself. Concentrations in
\f$[\text{threshold}, 0)\f$ are treated as floating-point noise (e.g. from MPI
partitioning or non-monotone advection/diffusion schemes) and clamped silently,
whereas concentrations \f$< \text{threshold}\f$ are clamped and reported via an
aggregated warning per speciation step. Must not be positive. Default: -1e-12
mol/kgw (about 0.002 ppb H2, close to laboratory detection limits).
