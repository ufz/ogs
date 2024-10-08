\copydoc MaterialPropertyLib::SaturationVanGenuchten

Either the `exponent` (which is the pressure exponent \f$m\f$ in the equations)
must be set and then the saturation exponent is \f$n = 1 / (1 - m)\f$.
Or both values the `pressure_exponent` and the `saturation_exponent` must be set
independently.
