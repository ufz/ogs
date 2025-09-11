# Collection of python boundary condition (BC) classes for OpenGeoSys

try:
    import ogs.callbacks as OpenGeoSys
except ModuleNotFoundError:
    import OpenGeoSys


s_a = 365.25 * 24 * 3600  # =31557600 seconds per year

# Choose parametrization
g = 9.81  # m/s²
nu = 0.4  # -
phi = 0.2  # -
rho_W = 1000.0  # kg/m³
rho_S = 3000.0  # kg/m³
u_max = -1.0  # m

Lx = 10  # m
Ly = 10  # m
T = 1e6  # s

# Nomenclature: BC Process_LocationQuantity_Component
# 					(THM)			(XYZ)

# Process	Dirichlet BC	Neumann BC (normal to boundary)
# T			temperature		heat flux
# H			pressure		hydraulic flux
# M			displacement	momentum flux (stress vector)


def ExternalDisplacement(x, t):
    return u_max * (t / T) * (x / Lx) ** 2


# Hydraulic BCs
# -------------
class BCH_SurfacePressure(OpenGeoSys.BoundaryCondition):
    def getDirichletBCValue(self, _t, coords, _node_id, _primary_vars):
        _x, _y, _z = coords

        value = 0.0

        return (True, value)


class BCH_VerticalPressure(OpenGeoSys.BoundaryCondition):
    def getDirichletBCValue(self, _t, coords, _node_id, _primary_vars):
        _x, y, _z = coords

        value = -rho_W * g * y

        return (True, value)


# Mechanics BCs
# -------------
class BCM_VerticalTraction_X(OpenGeoSys.BoundaryCondition):
    def getFlux(
        self, _t, coords, _primary_vars
    ):  # here Neumann BC: flux of linear momentum
        _x, y, _z = coords

        value = (
            (nu / (1 - nu) * ((phi - 1) * rho_W + (1 - phi) * rho_S) + rho_W) * g * y
        )
        jacob = [0.0, 0.0, 0.0]

        return (True, value, jacob)


class BCM_MonitoringTraction(OpenGeoSys.BoundaryCondition):
    def getFlux(self, t, coords, primary_vars):
        x, _y, _z = coords
        uy_actual = primary_vars[2]

        uy_target = ExternalDisplacement(x, t)
        uy_diff = uy_actual - uy_target

        if abs(uy_diff) > 1e-15:
            print("WARNING: uy_diff = ", uy_diff)

        value = 0.0
        jacob = [0.0, 0.0, 0.0]

        return (False, value, jacob)


class BCM_SurfaceDisplacement_Y(OpenGeoSys.BoundaryCondition):
    def getDirichletBCValue(self, t, coords, _node_id, _primary_vars):
        x, _y, _z = coords

        # prescribe displacement u_y
        value = ExternalDisplacement(x, t)

        return (True, value)


class BCM_BottomDisplacement_Y(OpenGeoSys.BoundaryCondition):
    def getDirichletBCValue(self, _t, _coords, _node_id, _primary_vars):
        # prescribe displacement u_y
        value = 0

        return (True, value)


class BCM_LateralDisplacement_X(OpenGeoSys.BoundaryCondition):
    def getDirichletBCValue(self, _t, _coords, _node_id, _primary_vars):
        # prescribe displacement u_x
        value = 0

        return (True, value)


# instantiate the BC objects used by OpenGeoSys
# ---------------------------------------------

AmbientPressure_py = BCH_SurfacePressure()
PressureProfile_py = BCH_VerticalPressure()

NormalStressProfile_py = BCM_VerticalTraction_X()
Traction4Monitoring_py = BCM_MonitoringTraction()

LateralDisplacement_py = BCM_LateralDisplacement_X()
SurfaceDisplacement_py = BCM_SurfaceDisplacement_Y()
BottomDisplacement_py = BCM_BottomDisplacement_Y()
