This flag enables two optimizations for the ComponentTransport process:

1. the global equation systems are assembled only once, in the first timestep,
   and reused subsequently
2. the non-linear solver will converge in a single iteration in each timestep

\attention
With this option enabled, **OGS will produce correct results only if the process equations are linear**.
OGS will not detect any non-linearities. It is the responsibility of the user to
ensure that the assembled equation systems are linear, indeed!
Furthermore, the material properties and parameters used in the process's
equations must **not carry any time dependence**.
