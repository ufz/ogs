###
# Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
# Distributed under a Modified BSD License.
# See accompanying file LICENSE.txt or
# http://www.opengeosys.org/project/license
###

import sys
from pathlib import Path

from pandas import read_csv

print(sys.version)
try:
    import ogs.callbacks as OpenGeoSys
except ModuleNotFoundError:
    import OpenGeoSys

df_server = read_csv(
    "initial.csv", delimiter=";", index_col=[0], dtype={"data_index": str}
)


def get_Tin(t):
    df_readfile = read_csv("readfile.txt", delimiter=";")
    time_list = df_readfile["time"].tolist()  # prepare the time adjustment
    t = min(
        time_list, key=lambda x: abs(x - t)
    )  # !!!makes simulation much slower!!! - time adjustment to nearest value in list in case time value isn't in the list (happens in Beier test)
    return [float(df_readfile.Tin[df_readfile.time == t])]


# OGS setting
# Dirichlet BCs
class BC(OpenGeoSys.BHENetwork):
    def initializeDataContainer(self):
        # initialize network and get data from the network
        # convert dataframe to column list
        t = 0  # 'initial time'
        data_col_1 = df_server[
            "Tin_val"
        ].tolist()  # 'Tin_val'; option: give array directly
        data_col_2 = df_server["Tout_val"].tolist()  # 'Tout_val'
        data_col_3 = df_server["Tout_node_id"].astype(int).tolist()  # 'Tout_node_id'
        data_col_4 = df_server["flowrate"].tolist()  # 'BHE flow rate'
        return (t, data_col_1, data_col_2, data_col_3, data_col_4)

    def serverCommunicationPreTimestep(self, t, _dt, Tin_val, Tout_val, flowrate):
        # TODO: Code for SimualtionX simulation
        # with t; only take the last results for each time point
        # TODO: say SimulationX the next time point from OGS

        Tin_val = get_Tin(t)
        with Path("T_in.txt").open(mode="a") as fd:
            fd.write(str(t) + str(Tin_val) + "\n")

        with Path("T_out.txt").open(mode="a") as fd:
            fd.write(str(t) + str(Tout_val) + "\n")

        with Path("flowrate.txt").open(mode="a") as fd:
            fd.write(str(t) + str(flowrate) + "\n")

        return (Tin_val, flowrate)

    def serverCommunicationPostTimestep(self, t, _dt, Tin_val, Tout_val, flowrate):
        Tin_val = [305]

        with Path("T_in.txt").open(mode="a") as fd:
            fd.write("post: " + str(t) + str(Tin_val) + "\n")

        with Path("T_out.txt").open(mode="a") as fd:
            fd.write("post: " + str(t) + str(Tout_val) + "\n")

        with Path("flowrate.txt").open(mode="a") as fd:
            fd.write("post: " + str(t) + str(flowrate) + "\n")


# main
bc_bhe = BC()
