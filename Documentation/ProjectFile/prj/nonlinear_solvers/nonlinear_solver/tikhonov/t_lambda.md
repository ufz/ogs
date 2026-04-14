Tikhonov lambda value that is added to the diagonal of the Jacobian matrix. Should be small and positive (typically in the range of 1e-20 to 1e-10).

A larger value provides stronger regularisation but may affect the accuracy of the solution. A smaller value provides weaker regularisation but may not be sufficient to improve convergence.

**Example:**

```xml
<lambda>1e-14</lambda>
```
