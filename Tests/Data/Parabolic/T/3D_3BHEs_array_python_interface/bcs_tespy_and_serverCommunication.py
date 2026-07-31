# SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
# SPDX-License-Identifier: BSD-3-Clause

import sys
from pathlib import Path

# TESPy imports CoolProp at module level, which segfaults inside OGS's
# embedded Python, so the workarounds have to be installed before TESPy is
# imported anywhere. See bhe_tespy_compat.py for the details.
import bhe_tespy_compat
import numpy as np
from pandas import read_csv

# The network model is read relative to the project file, not to the working
# directory: the ctest runs OGS in the build directory so that the log files
# written below do not end up in the source tree.
prj_dir = Path(ogs_prj_directory)  # noqa: F821

bhe_tespy_compat.install_workarounds()

from tespy.networks import Network  # noqa: E402

try:
    import ogs.callbacks as OpenGeoSys
except ModuleNotFoundError:
    import OpenGeoSys


print(sys.version)

# User setting ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# parameters
# refrigerant density
rho_f = 992.92  # kg/m3
# switch for special boundary conditions
# switch of the function for manually specified dynamic flowrate
switch_dyn_frate = "off"  # 'on','off'
# switch of the function for manually specified dynamic thermal demand
switch_dyn_demand = "on"  # 'on','off'
# consumer component name in the network model
consumer_name = "consumer"


# network status setting
def network_status(t):
    nw_status = "on"
    # month for closed network
    timerange_nw_off_month = []  # No month for closed network
    # t-1 to avoid the calculation problem at special time point,
    # e.g. t = 2592000.
    t_trans = int((t - 1) / 86400 / 30) + 1
    t_trans_month = t_trans
    if t_trans_month > 12:
        t_trans_month = t_trans - 12 * (int(t_trans / 12))
    if t_trans_month in timerange_nw_off_month:
        nw_status = "off"
    return nw_status


# dynamic consumer thermal load
def consumer_demand(t):  # dynamic thermal demand from consumer
    # time conversion
    t_trans = int((t - 1) / 86400 / 30) + 1
    if t_trans > 12:
        t_trans = t_trans - 12 * (int(t_trans / 12))
    # thermal demand in each month (assumed specific heat extraction rate*
    # length of BHE* number of BHE)
    month_demand = [
        -25 * 50 * 3,
        -25 * 50 * 3,
        -25 * 50 * 3,
        -25 * 50 * 3,
        -25 * 50 * 3,
        -25 * 50 * 3,
        -25 * 50 * 3,
        -25 * 50 * 3,
        -25 * 50 * 3,
        -25 * 50 * 3,
        -25 * 50 * 3,
        -25 * 50 * 3,
    ]
    return month_demand[t_trans - 1]


# dynamic hydraulic flow rate at the network inlet
def dyn_frate(t):
    # time conversion
    t_trans = int((t - 1) / 86400 / 30) + 1
    if t_trans > 12:
        t_trans = t_trans - 12 * (int(t_trans / 12))
    # flow rate in kg / s time curve in month
    month_frate = []
    return month_frate[t_trans - 1]


# End User setting +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


# create network dataframe
def create_dataframe():
    # return dataframe
    return read_csv(
        prj_dir / "pre" / "bhe_network.csv",
        delimiter=";",
        index_col=[0],
        dtype={"data_index": str},
    )


# TESPy calculation process
def get_tespy_results(t):
    # bhe network boundary conditions re parametrization
    # if network exist dynamic flowrate
    if switch_dyn_frate == "on":
        cur_frate = dyn_frate(t)
        localVars["inlet_name"].set_attr(m=cur_frate)
    # if network exist dynamic thermal demand
    if switch_dyn_demand == "on":
        # consumer thermal load:
        cur_month_demand = consumer_demand(t)
        consumer.set_attr(Q=cur_month_demand)
    # T_out re parametrization:
    for i in range(n_BHE):
        localVars["outlet_BHE" + str(i + 1)].set_attr(
            T=df.loc[data_index[i], "Tout_val"]
        )
    # solving network
    nw.solve(mode="design")
    # get Tin_val and flow rate
    for i in range(n_BHE):
        # get Tin_val
        df.loc[df.index[i], "Tin_val"] = (
            localVars["inlet_BHE" + str(i + 1)].get_attr("T").val
        )
        # get flowrate
        df.loc[df.index[i], "flowrate"] = (
            localVars["inlet_BHE" + str(i + 1)].get_attr("m").val / rho_f
        )
    return df["Tin_val"].tolist(), df["flowrate"].tolist()


# OGS setting
# Dirichlet BCs
class BC(OpenGeoSys.BHENetwork):
    def initializeDataContainer(self):
        # initialize network and get data from the network
        nw.solve(mode="design")
        get_tespy_results(0)
        # convert dataframe to column list
        t = 0  # 'initial time'
        data_col_1 = df["Tin_val"].tolist()  # 'Tin_val'
        data_col_2 = df["Tout_val"].tolist()  # 'Tout_val'
        data_col_3 = df["Tout_node_id"].astype(int).tolist()  # 'Tout_node_id'
        data_col_4 = df["flowrate"].tolist()  # 'BHE flow rate'
        return (t, data_col_1, data_col_2, data_col_3, data_col_4)

    def tespySolver(self, t, Tin_val, Tout_val):
        # network status:
        nw_status = network_status(t)
        # if network closed:
        if nw_status == "off":
            df.loc[:, "flowrate"] = 0
            cur_flowrate = df["flowrate"].tolist()
            return (True, True, Tout_val, cur_flowrate)

        # read Tout_val to dataframe
        for i in range(n_BHE):
            df.loc[df.index[i], "Tout_val"] = Tout_val[i]
        # TESPy solver
        cur_Tin_val, cur_flowrate = get_tespy_results(t)
        # check norm if network achieves the converge
        if_success = False
        pre_Tin_val = Tin_val
        norm = np.linalg.norm(abs(np.asarray(pre_Tin_val) - np.asarray(cur_Tin_val)))
        if norm < 10e-6:
            if_success = True
        # return to OGS
        return (True, if_success, cur_Tin_val, cur_flowrate)

    def serverCommunicationPreTimestep(self, t, _dt, Tin_val, Tout_val, flowrate):
        # TODO: Code for SimualtionX simulation
        # with t; only take the last results for each time point
        # TODO: say SimulationX the next time point from OGS

        with Path("T_in.txt").open(mode="a") as fd:
            fd.write(str(t) + str(Tin_val) + "\n")

        with Path("T_out.txt").open(mode="a") as fd:
            fd.write(str(t) + str(Tout_val) + "\n")

        with Path("flowrate.txt").open(mode="a") as fd:
            fd.write(str(t) + str(flowrate) + "\n")

        return (Tin_val, flowrate)

    def serverCommunicationPostTimestep(self, t, _dt, Tin_val, Tout_val, flowrate):
        Tin_val = [305, 305, 305]

        with Path("T_in.txt").open(mode="a") as fd:
            fd.write("post: " + str(t) + str(Tin_val) + "\n")

        with Path("T_out.txt").open(mode="a") as fd:
            fd.write("post: " + str(t) + str(Tout_val) + "\n")

        with Path("flowrate.txt").open(mode="a") as fd:
            fd.write("post: " + str(t) + str(flowrate) + "\n")

        return


# main
# initialize the tespy model of the bhe network
# loading the TESPy model
nw = Network.from_json(str(prj_dir / "pre" / "tespy_nw.json"))
# set if print the network iteration info
nw.iterinfo = False

# fails early and with a readable message if the label above was mistyped or
# the component was renamed in pre/3bhes.py
consumer = bhe_tespy_compat.get_component(nw, consumer_name)

# create bhe dataframe of the network system from bhe_network.csv
df = create_dataframe()
n_BHE = np.size(df.iloc[:, 0])

# create local variables of the components label and connections label in
# network
localVars = locals()
data_index = df.index.tolist()
for i in range(n_BHE):
    for c in nw.conns["object"]:
        # bhe inlet and outlet conns
        if c.target.label == data_index[i]:  # inlet conns of bhe
            localVars["inlet_BHE" + str(i + 1)] = c
        if c.source.label == data_index[i]:  # outlet conns of bhe
            localVars["outlet_BHE" + str(i + 1)] = c

# time depended flowrate
if switch_dyn_frate == "on":
    # the connection leaving the network's Source, i.e. the network inlet
    localVars["inlet_name"] = bhe_tespy_compat.find_network_inlet(nw)

# instantiate BC objects referenced in OpenGeoSys
bc_bhe = BC()
