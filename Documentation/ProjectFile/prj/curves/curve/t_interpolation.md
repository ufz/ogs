Interpolation between the supporting points of the curve. The default is `linear`.

- `linear`: piecewise linear interpolation. The derivative is not the exact derivative of the interpolant but an approximation: on the first and the last two intervals it is the slope of the respective interval, on the remaining intervals it blends the difference quotients of neighbouring supporting points linearly.
- `pchip`: shape-preserving piecewise cubic Hermite interpolation (`boost::math::interpolators::pchip`). The curve is continuously differentiable within the range of the supporting points, passes through all supporting points, and is monotone wherever the supporting values are monotone, i.e., it does not over- or undershoot the data. At least four supporting points are required.

Outside of the range of the supporting points the curve is constant for both interpolations. Hence, the derivative jumps to 0 at the first and last supporting points unless the respective end slope is 0.
