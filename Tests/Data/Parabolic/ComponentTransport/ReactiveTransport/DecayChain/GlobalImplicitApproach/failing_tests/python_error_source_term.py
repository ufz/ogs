try:
    import ogs.callbacks as OpenGeoSys
except ModuleNotFoundError:
    import OpenGeoSys


class ST(OpenGeoSys.SourceTerm):
    def getFlux(self, _t, _coords, _primary_vars):
        msg = "this exception is thrown on purpose"
        raise RuntimeError(msg)


st = ST()
