Thermal conductivity of the inlet pipe wall [W/(m*K)].
Accepts either a numeric value or a parameter name referencing a `<parameter>` defined in the project-level `<parameters>` block.
When a parameter name is given, the parameter is sampled at each BHE element's centroid to produce per-element wall conductivity values.
Thermal resistances are then computed independently for each element.
