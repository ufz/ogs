try:
    import ogs.callbacks as OpenGeoSys
except ModuleNotFoundError:
    import OpenGeoSys


class BC(OpenGeoSys.BoundaryCondition):
    def getDirichletBCValue(self, _t, _coords, _node_id, _primary_vars):
        msg = "this exception is thrown on purpose"
        raise RuntimeError(msg)


bc = BC()
