import OpenGeoSys



p_flux_in = 1e-2
p_0 = 1e5

# Source Terms
## Pressure
class BC_p_D(OpenGeoSys.BoundaryCondition):
    def getDirichletBCValue(self, t, coords, node_id, primary_vars):
        return (True, p_0)


class BC_p_N(OpenGeoSys.BoundaryCondition):
    def getFlux(self, t, coords, primary_vars):
        Jac = [0.0, 0.0, 0.0]
        return (True, p_flux_in, Jac)

bc_p_D = BC_p_D()
bc_p_N = BC_p_N()
