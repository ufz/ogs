try:
    import ogs.callbacks as OpenGeoSys
except ModuleNotFoundError:
    import OpenGeoSys
import chamber


class BC_X(OpenGeoSys.BoundaryCondition):
    def getDirichletBCValue(self, _t, coords, _node_id, _primary_vars):
        x, _y, _z = coords

        if x == 0.0 or x == chamber.R_piston:
            # no x displacement at outer boundary (x==R_piston) and center (x==0.0)
            return (True, 0.0)

        return (False, 0.0)


class BC_Y(OpenGeoSys.BoundaryCondition):
    def getDirichletBCValue(self, t, coords, _node_id, _primary_vars):
        _x, y, _z = coords

        if y == 0.04:
            # upper boundary -- prescribe displacement
            return (True, -t * chamber.h_chamber_0 / 20.0)

        return (False, 0.0)

    def getFlux(self, _t, coords, primary_vars):
        _x, y, _z = coords

        if y == 0.0:
            # lower boundary -- the chamber
            _ux, uy = primary_vars

            p = chamber.p_chamber(uy)
            Jac = [0.0, chamber.dp_chamber_du_y(uy)]
            return (True, p, Jac)

        return (False, 0.0, [0.0, 0.0])


# instantiate the BC objects used by OpenGeoSys
bc_x = BC_X()
bc_y = BC_Y()
