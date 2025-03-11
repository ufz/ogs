Parameter \f$\alpha\f$ for non-negative damping of the nonlinear solver.
The final damping \f$\delta\f$ is defined by the original damping \f$\delta_0\f$ with
 \f$\delta = \min_i \frac{\delta_0}{max(1.0, -\frac{\alpha_j \Delta x_i}{\cdot x_i}}\f$,
determined over all degrees of freedom $i$ and global components \f$j\f$.
A value of zero switches off the non-negative damping. A value of two, halves the step
if a given proposed step would cross the zero-line for a given component.
Higher values lead to a stronger damping.
