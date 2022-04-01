import OpenGeoSys


dirichlet_displacement_top = -0.05
dirichlet_displacement_0 = 0
dirichlet_pressure_0 = 0



# Boundary conditions
## Pressure
class BC_p_D_0(OpenGeoSys.BoundaryCondition):        
    def getDirichletBCValue(self, t, coords, node_id, primary_vars):
        return (True, dirichlet_pressure_0)
        
## Displacement
### Dirichlet
class BC_u_D_0(OpenGeoSys.BoundaryCondition):        
    def getDirichletBCValue(self, t, coords, node_id, primary_vars):
        return (True, dirichlet_displacement_0)
        
class BC_u_D_top(OpenGeoSys.BoundaryCondition):        
    def getDirichletBCValue(self, t, coords, node_id, primary_vars):
        return (True, dirichlet_displacement_top)


bc_u_D_0 = BC_u_D_0()
bc_u_D_top = BC_u_D_top()
bc_p_D_0 = BC_p_D_0()
