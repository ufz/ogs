Specific body forces applied to fluid.

It is usually used to apply gravitational forces.

The value is a single signed component along the vertical axis, in m/s^2, so
gravity is given as `-9.81`. The assembly projects it onto the well axis, that
is it multiplies it by the cosine of the inclination of the element measured
from the vertical.

Note that only this body force is projected onto the well axis. The drift-flux
closure the process solves the vapour void fraction with is a correlation for
vertical flow and receives no such correction, so a simulation on an inclined
mesh combines an inclination-corrected body force with an uncorrected buoyant
slip and is incorrect. The process is restricted to vertical wells; see the
documentation of `ProcessLib::WellboreSimulator::WellboreSimulatorProcess`.

A vanishing body force does not switch gravity off everywhere either: the drift
of the closure stays gravity-driven, so a gravity-free run is not internally
consistent.
