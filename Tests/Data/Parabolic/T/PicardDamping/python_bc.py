try:
    import ogs.callbacks as OpenGeoSys
except ModuleNotFoundError:
    import OpenGeoSys


class NonlinearRobinBC(OpenGeoSys.BoundaryCondition):
    def __init__(self):
        super().__init__()
        self.T_inf = 20.0
        self.h0 = 100.0

    def getFlux(self, _t, _coords, primary_vars):
        T = primary_vars[0]
        alpha = self.h0 / (1.0 + T * T)
        flux = alpha * (T - self.T_inf)
        dalpha_dT = -self.h0 * 2.0 * T / ((1.0 + T * T) ** 2)
        dFlux_dT = alpha + (T - self.T_inf) * dalpha_dT
        return (True, flux, [dFlux_dT])


nonlinear_robin_bc = NonlinearRobinBC()
