The borehole diameter.
Accepts either a numeric value (constant diameter) or a parameter name referencing a `<parameter>` defined in the project-level `<parameters>` block.
When a parameter name is given, the parameter is sampled at each BHE element's centroid to produce per-element diameter values (e.g., using a `Function` parameter with a step-function expression for telescoping boreholes).
Thermal resistances are then computed independently for each element.
