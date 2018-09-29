import OpenGeoSys


class BCPressure(OpenGeoSys.BoundaryCondition):
    def getDirichletBCValue(self, t, coords, node_id, primary_vars):
        x, y, z = coords

        if x == -150 and z == 75:
            # prescribe pressure of 0
            return (True, 0.0)

        # no Dirichlet BC
        return (False, 0.0)


class BCConcentration(OpenGeoSys.BoundaryCondition):
    def getDirichletBCValue(self, t, coords, node_id, primary_vars):
        x, y, z = coords

        if z == -75:
            # prescribe concentration of 0
            return (True, 0.0)

        if z == 75 and x >= 0:
            # prescribe concentration of 1
            return (True, 1.0)

        # no Dirichlet BC
        return (False, 0.0)


# instantiate the BC objects used by OpenGeoSys
bc_p = BCPressure()
bc_c = BCConcentration()
