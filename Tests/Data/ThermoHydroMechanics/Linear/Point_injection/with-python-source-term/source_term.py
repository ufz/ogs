try:
    import ogs.callbacks as OpenGeoSys
except ModuleNotFoundError:
    import OpenGeoSys


class SourceTerm(OpenGeoSys.SourceTerm):
    def getFlux(self, _t, coords, _primary_vars):
        x, y, z = coords
        value = 150
        Jac = [0.0] * 4
        return (value, Jac)


source_term = SourceTerm()
