A Jacobian assembler that assembles the Jacobian in two different ways, compares
the resulting local Jacobians and writes extensive logs in the form of a Python
script if the provided tolerances are exceeded.

Logging (and optionally program termination) is triggered only if both the
absolute and the relative tolerance are exceeded.
