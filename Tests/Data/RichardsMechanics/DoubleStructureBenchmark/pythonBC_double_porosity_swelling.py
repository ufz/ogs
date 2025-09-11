# Python boundary condition for a Dirichlet Jump with subsequent freeing of the boundary
try:
    import ogs.callbacks as OpenGeoSys
except ModuleNotFoundError:
    import OpenGeoSys


s_a = 365.25 * 24 * 3600  # =31557600 seconds per year

# Choose parametrization
t0 = 0.00  # s
tJ = 1.0e4  # s

p0 = 1e8  # Pa
pJ = -4e8  # Pa

# Nomenclature: BC Process_LocationQuantity_Component
# 					(THM)			(XYZ)

# Process	Dirichlet BC	Neumann BC (normal to boundary)
# H			pressure		hydraulic flux
# M			displacement	momentum flux (stress vector)


# Hydraulic BCs
# -------------
class BCH_PressureJump(OpenGeoSys.BoundaryCondition):
    def getDirichletBCValue(self, t, coords, _node_id, _primary_vars):
        _x, _y, _z = coords

        if t0 <= t <= tJ:
            value = p0 + (pJ - p0) * (t - t0) / (tJ - t0)
            return (True, value)
        return (False, 0)


# instantiate the BC objects used by OpenGeoSys
# ---------------------------------------------

InitialPressureJump_py = BCH_PressureJump()
