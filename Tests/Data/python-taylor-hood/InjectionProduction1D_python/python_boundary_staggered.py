import OpenGeoSys

dirichlet_displacement_0 = 0
neumann_displacement_overburden = -2.125e6
source_term_injection = 1.16e-4
source_term_production = -1.16e-4

# Source Terms
## Pressure
class Source_p_injection(OpenGeoSys.SourceTerm):
    def getFlux(self, t, coords, primary_vars):
        Jac = [0.0]
        return (source_term_injection, Jac)

class Source_p_production(OpenGeoSys.SourceTerm):
    def getFlux(self, t, coords, primary_vars):
        Jac = [0.0]
        return (source_term_production, Jac)

# Boundary conditions
## Pressure

## Displacement
### Dirichlet
class BC_u_D(OpenGeoSys.BoundaryCondition):        
    def getDirichletBCValue(self, t, coords, node_id, primary_vars):
        return (True, dirichlet_displacement_0)

class BC_u_N(OpenGeoSys.BoundaryCondition):
    def getFlux(self, t, coords, primary_vars):
        Jac = [0.0, 0.0]
        return (True, neumann_displacement_overburden, Jac)

injection = Source_p_injection()
production = Source_p_production()
bc_u_D = BC_u_D()
bc_u_N = BC_u_N()
