import OpenGeoSys
from chamber import *


class BC_X(OpenGeoSys.BoundaryCondition):
    def getDirichletBCValue(self, t, coords, node_id, primary_vars):
        x, y, z = coords

        if x == 0.0 or x == R_piston:
            # no x displacement at outer boundary (x==R_piston) and center (x==0.0)
            return (True, 0.0)

        return (False, 0.0)


class BC_Y(OpenGeoSys.BoundaryCondition):
    def getDirichletBCValue(self, t, coords, node_id, primary_vars):
        x, y, z = coords

        if y == 0.04:
            # upper boundary -- prescribe displacement
            return (True, -t * h_chamber_0 / 20.0)

        return (False, 0.0)


    def getFlux(self, t, coords, primary_vars):
        x, y, z = coords

        if y == 0.0:
            # lower boundary -- the chamber
            ux, uy = primary_vars

            p = p_chamber(uy)
            Jac = [ 0.0, dp_chamber_du_y(uy) ]
            return (True, p, Jac)

        return (False, 0.0, [ 0.0, 0.0 ])


# instantiate the BC objects used by OpenGeoSys
bc_x = BC_X()
bc_y = BC_Y()
