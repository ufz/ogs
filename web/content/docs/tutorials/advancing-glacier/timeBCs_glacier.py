import OpenGeoSys
import glacierclass as glc

L_dom = 120000  # m
L_max = 0.7 * L_dom  # m
H_max = 200  # m
x_0 = -0.5 * L_dom  # m
t_0 = 0.00  # s
t_1 = 1.0000  # s


class BC_Y(OpenGeoSys.BoundaryCondition):
    def __init__(self, L_dom, L_max, H_max, x_0, t_0, t_1):
        super(BC_Y, self).__init__()
        # instantiate the glacier member object
        self.glacier = glc.glacier(L_dom, L_max, H_max, x_0, t_0, t_1)
        self.glacier.printMaxLoads()
        self.glacier.plotEvolution()

    def getFlux(
        self, t, coords, primary_vars
    ):  # here Neumann BC: flux of linear momentum
        x, y, z = coords

        # print("x = ", x)
        # print("t = ", t)
        # print("s = ", self.glacier.normalstress_glacier(x,t))
        value = self.glacier.normalstress(x, t)
        derivative = [0.0, 0.0]

        return (True, value, derivative)


# instantiate the BC objects used by OpenGeoSys
bc_y = BC_Y(L_dom, L_max, H_max, x_0, t_0, t_1)
