import numpy as np
from AfuncW import NGP, CIC, TSC, PSQ, W
import random


class grid:
    def __init__(self, order):
        self.order = order
        if order == 0:
            self.func = NGP
        elif order == 1:
            self.func = CIC
        elif order == 2:
            self.func = TSC
        elif order == 3:
            self.func = PSQ
        else:
            print("There are just assignment functions up to order 3")
            print("implemented yet /n")

    def create_grid(self):
        """creates grid appropiate to the order of the assignement function

        returns meshgrids for each direction
        """
        order = self.order
        # number of cells, that a-func fits in
        numCells = order + 3 + 1  # +1 for recursive current calculation [i-1]
        # creates the grids axis
        axis = np.arange(0, numCells, 1, dtype=float)
        # creates grid, meshgrid for each direction
        grid_x, grid_y, grid_z = np. meshgrid(axis, axis, axis)

        return grid_x, grid_y, grid_z

    def ran_partikle_step(self):
        """Randomly selects start and end point of a macro particle.
        Randomly selects whether the endpoint is located in the initial cell
        or not;
        """
        order = self.order
        num_cells = order + 3 + 1

        # choose Startcell in the middle of the grid
        start_x_zelle = int((num_cells-1) / 2)
        start_y_zelle = int((num_cells-1) / 2)
        start_z_zelle = int((num_cells-1) / 2)

        # random
        start_x = start_x_zelle + random.random()
        start_y = start_y_zelle + random.random()
        start_z = start_z_zelle + random.random()

        start_coord = np.array([start_x, start_y, start_z])

        # choose endcoordinates
        end_x = start_x + random.uniform(-1, 1)
        end_y = start_y + random.uniform(-1, 1)
        end_z = start_z + random.uniform(-1, 1)

        end_coord = np.array([end_x, end_y, end_z])

        return start_coord, end_coord

    def partikle_step(self, pos1, pos2, pos_off1, pos_off2):
        """Performs a deterministic particle step for a particle in the
        minimal grid. It takes the coordinate tuples pos1 and pos2, which
        describe the particle's position within the cell, as well as
        pos_off1/2 as tuples for the node coordinates.
        """
        order = self.order
        num_cells = order + 3 + 1

        # Choose starting cell at the center of the grid
        start_x_cell = int((num_cells - 1) / 2)
        start_y_cell = int((num_cells - 1) / 2)
        start_z_cell = int((num_cells - 1) / 2)

        # Coordinates in the middle cell + position in the cell
        start_x = start_x_cell + pos1[0]
        start_y = start_y_cell + pos1[1]
        start_z = start_z_cell + pos1[2]

        start_coord = np.array([start_x, start_y, start_z])

        # Calculate end coordinates - starting from starting cell, add position
        # and any cell crossings pos_off1 - pos_off2
        end_x = start_x_cell + pos2[0] - pos_off1[0] + pos_off2[0]
        end_y = start_y_cell + pos2[1] - pos_off1[1] + pos_off2[1]
        end_z = start_z_cell + pos2[2] - pos_off1[2] + pos_off2[2]

        end_coord = np.array([end_x, end_y, end_z])

        return start_coord, end_coord

    def current_deposition_field(self, start_coord, end_coord, grid_x,
                                 grid_y, grid_z):
        """Calculates the current deposition vector component-wise for z, y, x.
        Returns the components of the vector in the x, y, z direction on the
        grid. Takes the tuple of start and end coordinates (x, y, z) as well
        as the zero grids for each coordinate.
        """
        order = self.order
        func = self.func
        num_cells = order + 3 + 1

        for z in range(num_cells):
            for y in range(num_cells):
                for x in range(num_cells):
                    # Calculate the si components for the W calculation

                    # Old
                    s1 = func(start_coord[0] - x)
                    s2 = func(start_coord[1] - y)
                    s3 = func(start_coord[2] - z)
                    # new
                    s4 = func(end_coord[0] - x)
                    s5 = func(end_coord[1] - y)
                    s6 = func(end_coord[2] - z)

                    # calculate W and assign it to the grid node
                    grid_x[z][y][x] = W(s1, s2, s3, s4, s5, s6)
                    grid_y[z][y][x] = W(s2, s1, s3, s5, s4, s6)
                    grid_z[z][y][x] = W(s3, s2, s1, s6, s5, s4)

        return grid_x, grid_y, grid_z

    def current_density_field(self, W_grid_x, W_grid_y, W_grid_z, start_coord,
                              end_coord, Q, cell_depth, cell_width,
                              cell_height, dt, unit_time, unit_charge,
                              unit_length):
        """Calculates the current density in the x, y, and z directions on the
        grid using the current deposition vector(field) W_grid_x,y,z. Charge Q
        is needed as input. dx, dy, dz, and dt are set to 1 by default.
        """

        order = self.order
        num_cells = order + 3 + 1
        c = 299792458.0

        # create null grid to write current density to
        # due to periodicity, the grid has the same shape -> i = 1/2 is grid[0]
        null_axis = np.zeros(num_cells, dtype=float)
        j_grid_x, j_grid_y, j_grid_z = np.meshgrid(null_axis, null_axis,
                                                   null_axis)

        # determine how many grid points need to be included in calculation
        # aka half-width of assignment function
        # coordinate of coordinate to be calculated
        add_active = (order + 1) * 1/2

        # determine start x
        start_i = int(min(start_coord[0], end_coord[0]) - add_active)
        start_j = int(min(start_coord[1], end_coord[1]) - add_active)
        start_k = int(min(start_coord[2], end_coord[2]) - add_active)

        # determine end x
        end_i = int(max(start_coord[0], end_coord[0]) + add_active) + 1
        end_j = int(max(start_coord[1], end_coord[1]) + add_active) + 1
        end_k = int(max(start_coord[2], end_coord[2]) + add_active) + 1

        x_factor = (-Q * unit_charge) / (c * dt * unit_time) * (cell_width *
                                                                unit_length)

        # calculate the norm of the assignment function
        H = cell_width * cell_height * cell_depth * unit_length**3
        norm = (1/H)**order

        for k in range(start_k, end_k, 1):
            for j in range(start_j, end_j, 1):
                # -1 at end_i to ensure that calculation stops at last W_x
                for i in range(start_i, end_i - 1, 1):
                    j_grid_x[k][j][i] = (x_factor * norm * W_grid_x[k][j][i] +
                                         j_grid_x[k][j][i-1])

        y_factor = (-Q * unit_charge) / (c * dt * unit_time) * (cell_height
                                                                * unit_length)
        for k in range(start_k, end_k, 1):
            # -1 at end_j to ensure that calculation stops at last W_y
            for j in range(start_j, end_j - 1, 1):
                for i in range(start_i, end_i, 1):
                    j_grid_y[k][j][i] = (y_factor * norm * W_grid_y[k][j][i]
                                         + j_grid_y[k][j-1][i])

        z_factor = (-Q * unit_charge) / (c * dt * unit_time) * (cell_depth *
                                                                unit_length)
        # -1 at end_k to ensure that calculation stops at last W_z
        for k in range(start_k, end_k - 1, 1):
            for j in range(start_j, end_j, 1):
                for i in range(start_i, end_i, 1):
                    j_grid_z[k][j][i] = (z_factor * norm * W_grid_z[k][j][i] +
                                         j_grid_y[k-1][j][i])

        return j_grid_x, j_grid_y, j_grid_z
