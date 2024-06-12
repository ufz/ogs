try:
    import ogs.callbacks as OpenGeoSys
except ModuleNotFoundError:
    import OpenGeoSys


class BC(OpenGeoSys.BoundaryCondition):
    def getFlux(self, _t, _coords, _primary_vars):
        msg = "this exception is thrown on purpose"
        raise RuntimeError(msg)


bc = BC()
