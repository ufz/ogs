from __future__ import print_function

import OpenGeoSys


# enable/disable debug output
if True:
    def debug(*args):
        print(*args)
else:
    def debug(*args):
        pass


SPHERE_RADIUS = 1.0
START_TIME = 0.0


class HertzContactBC(OpenGeoSys.BoundaryCondition):
    def __init__(self):
        super(HertzContactBC, self).__init__()

        self._first_node = None         # ID of the first node of this BC's geometry
        self._t_old = START_TIME - 1.0  # time of previous invocation of this BC


    def _init_timestep(self):
        """Initializes the internal state at the beginning of a new timestep."""

        self._a_range = [ 0.0, SPHERE_RADIUS ]  # range of possible contact radii
        self._a_est = SPHERE_RADIUS             # estimated contact radius

        self._max_x_with_y_excess = -1.0        # maximum x value where a node is above the contact line


    def _init_iteration(self):
        """Initializes the internal state at the beginning of a new nonlinear solver iteration."""

        # variant of a bisection algorithm
        if self._max_x_with_y_excess >= 0.0:
            # there were nodes above the contact area in the previous
            # iteration => contact area has to grow.
            if self._max_x_with_y_excess < self._a_est:
                # This is an ad-hoc variation of the bisection algorithm for
                # finding the radius of the contact area. The variation is
                # necessary, because we obtain the value of
                # self._max_x_with_y_excess only after we already applied the
                # BCs for the current nonlinear solver iteration, i.e., the
                # applied BCs and the information based on which the BCs are
                # applied are off by one iteration.
                self._a_range[1] = self._max_x_with_y_excess
                assert self._a_range[0] <= self._a_range[1]
            else:
                self._a_range[0] = self._a_est
        else:
            # contact area has to shrink
            self._a_range[1] = self._a_est

        self._a_est = 0.5 * sum(self._a_range)
        print("BC: a_est={:.4f}, ({:.4f}, {:.4f})".format(self._a_est, self._a_range[0], self._a_range[1]))

        self._max_x_with_y_excess = -1.0


    def getDirichletBCValue(self, t, coords, node_id, primary_vars):
        if self._t_old < t:
            self._t_old = t
            self._init_timestep()

        # detect nonlinear iteration, assumes that nodes are always processed in
        # the same order
        if self._first_node is None:
            self._first_node = node_id

        if self._first_node == node_id:
            self._init_iteration()

        x, y, z = coords
        ux, uy = primary_vars

        try:
            # check that we are at the outer boundary
            assert abs(x**2 + y**2 + z**2 - 1.0) < 1e-15
        except:
            print("assert abs(x**2 + y**2 + z**2 - 1.0) < 1e-15",
                    x, y, z, x**2 + y**2 + z**2 - 1)
            raise

        y_deformed = y + uy
        y_top = HertzContactBC._get_y_top(t)

        res = (False, 0.0)

        if y_deformed >= y_top:
            self._max_x_with_y_excess = max(self._max_x_with_y_excess, x)

            if x <= self._a_est:
                res = (True, y_top - y)

        return res


    @staticmethod
    def _get_y_top(t):
        """Returns the y-position of the contact area depending on the load step t."""

        return 1.0 - 0.005 * t


# instantiate the BC object used by OpenGeoSys
bc = HertzContactBC()
