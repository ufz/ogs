
import OpenGeoSys
from math import pi, sin, cos, sinh, cosh

a = 2.0*pi/3.0

# analytical solution used to set the Dirichlet BCs
def solution(x, y):
    return sin(a*x) * sinh(a*y)

# gradient of the analytical solution used to set the Neumann BCs
def grad_solution(x, y):
    return a * cos(a*x) * sinh(a*y), \
            a * sin(a*x) * cosh(a*y)

# Dirichlet BCs
class BCTop(OpenGeoSys.BoundaryCondition):
    def getDirichletBCValue(self, t, coords, node_id, primary_vars):
        x, y, z = coords
        assert y == 1.0 and z == 0.0
        value = solution(x, y)
        return (True, value)

class BCLeft(OpenGeoSys.BoundaryCondition):
    def getDirichletBCValue(self, t, coords, node_id, primary_vars):
        x, y, z = coords
        assert x == 0.0 and z == 0.0
        value = solution(x, y)
        return (True, value)

class BCBottom(OpenGeoSys.BoundaryCondition):
    def getDirichletBCValue(self, t, coords, node_id, primary_vars):
        x, y, z = coords
        assert y == 0.0 and z == 0.0
        value = solution(x, y)
        return (True, value)

# Neumann BC
class BCRight(OpenGeoSys.BoundaryCondition):
    def getFlux(self, t, coords, primary_vars):
        x, y, z = coords
        assert x == 1.0 and z == 0.0
        value = grad_solution(x, y)[0]
        Jac = [ 0.0 ]  # value does not depend on primary variable
        return (True, value, Jac)


# instantiate BC objects referenced in OpenGeoSys' prj file
bc_top = BCTop()
bc_right = BCRight()
bc_bottom = BCBottom()
bc_left = BCLeft()
