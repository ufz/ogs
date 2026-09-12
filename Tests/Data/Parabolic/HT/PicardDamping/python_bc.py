try:
    import ogs.callbacks as OpenGeoSys
except ModuleNotFoundError:
    import OpenGeoSys


class NonlinearRobinBC(OpenGeoSys.BoundaryCondition):
    """Same destabilizing mechanism as Tests/Data/Parabolic/T/PicardDamping's
    python_bc.py (alpha = h0/(1+T^2), which weakens for large |T|, is exactly
    what forces plain Picard to oscillate there without damping/acceleration),
    shifted by T_ref so it pivots at the same point relative to Kelvin-scale
    temperatures as the original test's did relative to its 0-based ones.
    """

    def __init__(self):
        super().__init__()
        self.T_ref = 273.0
        self.T_inf = 20.0
        self.h0 = 100.0

    def getFlux(self, _t, _coords, primary_vars):
        # Monolithic HT: primary_vars is [T, p], and the Jacobian contribution
        # of this BC must carry a derivative entry for every primary variable
        # at the node, even though the flux itself does not depend on p.
        T = primary_vars[0] - self.T_ref
        alpha = self.h0 / (1.0 + T * T)
        flux = alpha * (T - self.T_inf)
        dalpha_dT = -self.h0 * 2.0 * T / ((1.0 + T * T) ** 2)
        dFlux_dT = alpha + (T - self.T_inf) * dalpha_dT
        return (True, flux, [dFlux_dT, 0.0])


nonlinear_robin_bc = NonlinearRobinBC()
