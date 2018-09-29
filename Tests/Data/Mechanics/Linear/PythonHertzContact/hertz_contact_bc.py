from __future__ import print_function

import OpenGeoSys


SPHERE_RADIUS = 1.0
START_TIME = 0.0


class HertzContactBC(OpenGeoSys.BoundaryCondition):
    def __init__(self):
        super(HertzContactBC, self).__init__()

        self._first_node = None         # ID of the first node of this BC's geometry
        self._t_old = START_TIME - 1.0  # time of previous invocation of this BC

        self._boundary_x_coords = []    # the x coordinates of all boundary nodes


    def _init_timestep(self):
        """Initializes the internal state at the beginning of a new timestep."""

        self._a_range = [ 0.0, SPHERE_RADIUS ]  # range of possible contact radii
        self._a_est = SPHERE_RADIUS             # estimated contact radius

        self._max_x_with_y_excess = -1.0        # maximum x value where a node is above the contact line

        self._tendency_is_up = True             # whether the contact area is supposed to grow in the current iteration
        self._a_curr = 0.0                      # the radius of the contact area in this iteration (continuously updated)
        self._a_prev = 0.0                      # the radius of the contact area that was prescripted in the previous iteration

        self._iteration = 0                     # the nonlinear solver iteration number


    def _init_iteration(self):
        """Initializes the internal state at the beginning of a new nonlinear solver iteration."""

        # variant of a bisection algorithm
        if self._max_x_with_y_excess >= 0.0:
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
                self._tendency_is_up = False
            else:
                # there were nodes above the contact area in the previous
                # iteration => contact area has to grow.
                self._a_range[0] = self._a_est
                self._tendency_is_up = True
        else:
            # contact area has to shrink
            self._a_range[1] = self._a_est
            self._tendency_is_up = False

        self._a_est = 0.5 * sum(self._a_range)

        self._max_x_with_y_excess = -1.0
        self._a_prev = self._a_curr
        self._a_curr = -1.0

        self._iteration += 1

        print("BC: a_est={:.4f}, a_prev={:.4f} ({:.4f}, {:.4f})".format(
            self._a_est, self._a_prev, self._a_range[0], self._a_range[1]))


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
            if self._iteration == 2:
                self._boundary_x_coords.sort()

        x, y, z = coords
        ux, uy = primary_vars

        if not self._boundary_x_coords_are_initialized(t):
            self._boundary_x_coords.append(x)

        try:
            # check that we are at the outer boundary
            assert abs(x**2 + y**2 + z**2 - SPHERE_RADIUS**2) < 1e-15
        except:
            print("assert abs(x**2 + y**2 + z**2 - 1.0) < 1e-15",
                    x, y, z, x**2 + y**2 + z**2 - SPHERE_RADIUS**2)
            raise

        y_deformed = y + uy
        y_top = HertzContactBC._get_y_top(t)

        res = (False, 0.0)

        if y_deformed >= y_top:
            self._max_x_with_y_excess = max(self._max_x_with_y_excess, x)

            if x <= self._a_est:
                res = (True, y_top - y)
                self._a_curr = max(x, self._a_curr)

            elif self._tendency_is_up:
                # This branch catches some corner-cases where Dirichlet BCs are
                # set on too few nodes.
                #
                # The existence of these corner cases and the fact that the
                # bisection algorithm for finding the correct contact radius is
                # rather ad-hoc show that in general for similar problems it
                # will be necessary to have more flexible convergence criteria
                # or nonlinear solvers.  In particular said similar problems are
                # those where the surfaces where Dirichlet or Neumann BCs are
                # set change for changing primary variables.

                if x <= self._a_prev:
                    assert False  # this case shouldn't happen
                    res = (True, y_top - y)
                elif self._boundary_x_coords_are_initialized(t):
                    idx = self._boundary_x_coords.index(x)
                    if idx != 0 and self._boundary_x_coords[idx-1] == self._a_prev:
                        res = (True, y_top - y)


        return res


    @staticmethod
    def _get_y_top(t):
        """Returns the y-position of the contact area depending on the load step t."""

        return 1.0 - 0.005 * t


    def _boundary_x_coords_are_initialized(self, t):
        return self._iteration >= 2


# instantiate the BC object used by OpenGeoSys
bc = HertzContactBC()
