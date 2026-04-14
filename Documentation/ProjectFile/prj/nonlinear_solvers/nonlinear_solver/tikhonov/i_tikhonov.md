Defines the configuration for Tikhonov regularization in the Newton nonlinear solver.

The Tikhonov regularization adds a value to the diagonal of the Jacobian matrix to improve convergence behaviour for ill-conditioned systems.

**Example:**

```xml
<tikhonov>
    <lambda>1e-14</lambda>
    <starting_iteration>11</starting_iteration>
</tikhonov>
```

The Tikhonov regularization is only applied starting from the iteration specified by `starting_iteration` to avoid interfering with the initial convergence phase.
