This parameter introduces a damping reduction \f$\lambda\f$ (\f$0 < \lambda < \f$ `inf`) by specifying the iteration number \f$i\f$ that it takes to reduce the damping \f$\delta\f$ to 1:
\f$ \delta + (1 - \delta) * \f$ `std::clamp(` \f$i / \lambda \f$, 0.0, 1.0`)`.
