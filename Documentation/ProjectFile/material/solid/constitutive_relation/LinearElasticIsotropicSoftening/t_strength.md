strength parameter governing the material's stress and stiffness.

Provide it as a time dependent parameter which starts as 1.0 (meaning
100% material strength) and ends at the target strength (e.g. 0.1 for
10% remaining material strength) to ensure a gradual softening.
The value corresponds to the ratio between the current stress trace and the
initial stress trace.
If the strength reaches 0, convergence might not be possible anymore. A close,
but non-zero target value might be more favorable to reach an almost
zero-strength state.
