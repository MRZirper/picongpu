import openpmd_api as io
import numpy as np
from grid_class import grid
import sys


def get_params(path):
    """creates a series of the data found at path.
    path should have the form »data_%T« to provide access to all iterations.

    returns: "./simData_%T.h5"
        series
        order... order of assignement function
        parameter list... list of parameters of the simulation (units etc.)"""
    # load series
    series = io.Series(path, io.Access.read_only)

    if len(series.iterations) == 1:
        print("There is just 1 iteration in the series")
        print("make sure, there are at least two")
        longer_series = False

    else:
        longer_series = True

    # read order assignment function
    i = series.iterations[0]
    electrons = i.particles["e"]
    order = int(electrons.get_attribute("particleShape"))

    # read parameters of the simulation
    cell_depth = i.get_attribute("cell_depth")
    cell_height = i.get_attribute("cell_height")
    cell_width = i.get_attribute("cell_width")
    dt = i.get_attribute("dt")
    unit_time = i.get_attribute("unit_time")
    unit_length = i.get_attribute("unit_length")
    unit_charge = i.get_attribute("unit_charge")

    return longer_series, series, order, [
        cell_depth, cell_height, cell_width, dt, unit_time, unit_charge,
        unit_length
        ]


def read_series(series):
    """reads the position values »position_offset« (grid_knot) and »position«
    (float in the intervall [0,1), relativ position in the cell) and the
    current density provided by the simulation as well the charge of the
    particle.

    returns:
        Js... current density list (contains array of each component for each
                                    iteration)
        poss... array of all relative in-Cell positions
        pos_offs... array of all knot-positions
        charge.. charge of particle
        """
    num_iter = len(series.iterations)

    # create arrays for saving the data
    Js = []
    poss = np.zeros((num_iter, 3))
    pos_offs = np.zeros((num_iter, 3))

    for iteration in range(0, num_iter, 1):
        i = series.iterations[iteration]  # current position

        # current density field
        J_x = i.meshes["J"]["x"].load_chunk()
        J_y = i.meshes["J"]["y"].load_chunk()
        J_z = i.meshes["J"]["z"].load_chunk()

        # charge
        charge = i.particles["e"]["charge"][io.Mesh_Record_Component.SCALAR]

        # node coordinates
        pos2_x = i.particles["e"]["position"]["x"].load_chunk()
        pos2_y = i.particles["e"]["position"]["y"].load_chunk()
        pos2_z = i.particles["e"]["position"]["z"].load_chunk()

        # InCell coordinates
        pos_off2_x = i.particles["e"]["positionOffset"]["x"].load_chunk()
        pos_off2_y = i.particles["e"]["positionOffset"]["y"].load_chunk()
        pos_off2_z = i.particles["e"]["positionOffset"]["z"].load_chunk()

        # Spühlen nicht vergessen!
        series.flush()

        # write coordinatetupel
        pos2 = np.array([*pos2_x, *pos2_y, *pos2_z])
        pos_off = np.array([*pos_off2_x, *pos_off2_y, *pos_off2_z])
        # write current density tupel
        J = np.array([J_x, J_y, J_z])

        Js.append(J)
        poss[iteration] = pos2
        pos_offs[iteration] = pos_off

        J_x_unit = i.meshes["J"]["x"]
        J_convert = J_x_unit.unit_SI

    return Js, poss, pos_offs, charge, J_convert


def compare(j_grid_x, j_grid_y, j_grid_z, J, J_convert):
    """compares the current desity provided by the simulation with the one
    calculated by the algorithm

    if both fields are identical, return 0
    if they are different in at least one component return 42"""

    c = 299792458.0  # LICHTGESCHWINDIGKEIT!!!""""
    # put all current desity components (unequal to 0) into a 1D Array
    x_sim_snake = J[0][J[0] != 0]
    x_algo_snake = j_grid_x[j_grid_x != 0]

    y_sim_snake = J[1][J[1] != 0]
    y_algo_snake = j_grid_y[j_grid_y != 0]

    z_sim_snake = J[2][J[2] != 0]
    z_algo_snake = j_grid_z[j_grid_z != 0]

    # compare the arrays values by deviding; -1 for np.where operation
    x_compare = x_sim_snake / x_algo_snake * J_convert / c - 1
    y_compare = y_sim_snake / y_algo_snake * J_convert / c - 1
    z_compare = z_sim_snake / z_algo_snake * J_convert / c - 1

    print("(algorithmarray / simulationsarray) ([converted])")
    print(x_compare + 1)
    print(65 * "=")

    epsilon = 1e-6
    x_compare = np.where(abs(x_compare) < epsilon, True, False)
    y_compare = np.where(abs(y_compare) < epsilon, True, False)
    z_compare = np.where(abs(z_compare) < epsilon, True, False)

    print("return of the boolean Array")
    print(x_compare)
    print(65 * "=")

    if np.all(x_compare) and np.all(y_compare) and np.all(z_compare):
        print("simulation and algorithm coincide")
        return 0
    else:
        print("no consensus between simulation and algorithm")
        print("please check the simulation")
        return 42


def main():

    path = sys.argv[1]  # Path for series, given as console argument
    # read parameters of series
    longer_series, series, order, params = get_params(path)

    if not longer_series:
        return 42

    # read Current deseties of series, charge, SI-convertion-factor, positions
    Js, poss, pos_offs, charge, J_convert = read_series(series)

    # create grid_object for comparison
    compare_grid = grid(order)

    # just for the first step (iteration 0&1)
    Q = charge[0]
    J = Js[1]
    pos1 = poss[0]
    pos2 = poss[1]
    pos_off1, pos_off2 = pos_offs[0], pos_offs[1]

    grid_x, grid_y, grid_z = compare_grid.create_grid()
    start_koord, end_koord = compare_grid.partikle_step(pos1, pos2,
                                                        pos_off1, pos_off2)
    # calculations via algorithm
    W_grid_x, W_grid_y, W_grid_z = (
        compare_grid.current_deposition_field(start_koord, end_koord, grid_x,
                                              grid_y, grid_z))

    j_grid_x, j_grid_y, j_grid_z = compare_grid.current_density_field(
        W_grid_x, W_grid_y, W_grid_z, start_koord, end_koord, Q, *params)

    # comparison
    compare_result = compare(j_grid_x, j_grid_y, j_grid_z, J, J_convert)

    return compare_result


if __name__ == "__main__":
    main()
