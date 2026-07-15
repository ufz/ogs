The iteration number from which Tikhonov regularisation starts to be applied. The default value is `0`, meaning the regularisation is applied from the first iteration.

Should be chosen such that the regularisation does not interfere with the initial convergence phase. Often set to a value where the solver shows poor convergence (e.g., 10-15 iterations).

**Example:**

```xml
<starting_iteration>11</starting_iteration>
```
