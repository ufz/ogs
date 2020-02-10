###
# Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
# Distributed under a Modified BSD License.
# See accompanying file LICENSE.txt or
# http://www.opengeosys.org/project/license
###

import sys
print(sys.version)
import os
import numpy as np
from pandas import read_csv
import OpenGeoSys
from tespy.networks import load_network

# User setting +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# parameters
# refrigerant parameters
refrig_density = 992.92  # kg/m3
# switch for special boundary conditions
# 'on','off', switch of the function for dynamic thermal demand from consumer
switch_dyn_demand = 'on'
# 'on','off', switch of the function for dynamic flowrate in BHE
switch_dyn_frate = 'off'


# network status setting
def network_status(t):
    nw_status = 'on'
    # month for closed network
    timerange_nw_off_month = [-9999]  # No month for closed network
    # t-1 to avoid the calculation problem at special time point,
    # e.g. t = 2592000.
    t_trans = int((t - 1) / 86400 / 30) + 1
    t_trans_month = t_trans
    if t_trans_month > 12:
        t_trans_month = t_trans - 12 * (int(t_trans / 12))
    if t_trans_month in timerange_nw_off_month:
        nw_status = 'off'
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
        -25 * 50 * 3, -25 * 50 * 3, -25 * 50 * 3, -25 * 50 * 3, -25 * 50 * 3,
        -25 * 50 * 3, -25 * 50 * 3, -25 * 50 * 3, -25 * 50 * 3, -25 * 50 * 3,
        -25 * 50 * 3, -25 * 50 * 3
    ]
    return month_demand[t_trans - 1]


# dynamic hydraulic flow rate
def dyn_frate(t):  # dynamic flowrate in BHE
    # time conversion
    t_trans = int((t - 1) / 86400 / 30) + 1
    if t_trans > 12:
        t_trans = t_trans - 12 * (int(t_trans / 12))
    # flow rate in kg / s time curve in month
    month_frate = [-9999]
    return month_frate[t_trans - 1]


# End User setting+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


# create network dataframe
def create_dataframe():
    # return dataframe
    df_nw = read_csv('./pre/bhe_network.csv',
                        delimiter=';',
                        index_col=[0],
                        dtype={'data_index': str})
    return (df_nw)


# TESPy hydraulic calculation process
def get_hydraulics(t):
    # if network exist dynamic flowrate
    if switch_dyn_frate == 'on':
        cur_frate = dyn_frate(t)
        localVars['inlet_name'].set_attr(m=cur_frate)
    # solve imported network
    nw.solve(mode='design')
    # get flowrate #m ^ 3 / s
    for i in range(n_BHE):
        for c in nw.conns.index:
            if c.t.label == data_index[i]:  # t:inlet comp, s:outlet comp
                df.loc[df.index[i], 'flowrate'] = c.get_attr('m').val_SI / refrig_density
    return df


# TESPy Thermal calculation process
def get_thermal(t):
    # bhe network thermal re parametrization
    if switch_dyn_demand == 'on':
        # consumer thermal load:
        cur_month_demand = consumer_demand(t)
        # print('cur_month_demand', cur_month_demand)
        nw.busses[bus_name].set_attr(P=cur_month_demand)
    # T_out:
    for i in range(n_BHE):
        localVars['outlet_BHE' + str(i + 1)].set_attr(T=df.loc[data_index[i],
                                                               'Tout_val'])
    # print('Tout=', df.loc[data_index[i], 'Tout_val'])
    # solving network
    nw.solve(mode='design')
    # get Tin_val
    for i in range(n_BHE):
        df.loc[df.index[i],
               'Tin_val'] = localVars['inlet_BHE' +
                                      str(i + 1)].get_attr('T').val_SI
    # print('Tin=', df.loc[df.index[i], 'Tin_val'])
    return df['Tin_val'].tolist()


# OGS setting
# Dirichlet BCs
class BC(OpenGeoSys.BHENetwork):
    def initializeDataContainer(self):
        # convert dataframe to column list
        t = 0  # 'initial time'
        data_col_1 = df['Tin_val'].tolist()  # 'Tin_val'
        data_col_2 = df['Tout_val'].tolist()  # 'Tout_val'
        data_col_3 = df['Tout_node_id'].astype(int).tolist()  # 'Tout_node_id'
        get_hydraulics(0)
        data_col_4 = df['flowrate'].tolist()  # 'BHE flow rate'
        return (t, data_col_1, data_col_2, data_col_3, data_col_4)

    def tespyThermalSolver(self, t, Tin_val, Tout_val):
        # network status:
        nw_status = network_status(t)
        # if network closed:
        #     print('nw_status = ', nw_status)
        if nw_status == 'off':
            return (True, True, Tout_val)
        else:
            # read Tout_val to dataframe
            for i in range(n_BHE):
                df.loc[df.index[i], 'Tout_val'] = Tout_val[i]
            # TESPy solver
            cur_cal_Tin_val = get_thermal(t)
            # check norm if network achieves the converge
            if_success = False
            pre_cal_Tin_val = Tin_val
            norm = np.linalg.norm(
                abs(np.asarray(pre_cal_Tin_val) - np.asarray(cur_cal_Tin_val)))
            if norm < 10e-6:
                if_success = True
            # return to OGS
            return (True, if_success, cur_cal_Tin_val)

    def tespyHydroSolver(self, t):
        if_dyn_frate = False
        data_flowrate = df['flowrate'].tolist()
        if switch_dyn_frate == 'on':
            if_dyn_frate = True
            # network status:
            nw_status = network_status(t)
            if nw_status == 'off':
                for i in range(n_BHE):
                    df.loc[df.index[i], 'flowrate'] = 0
                data_flowrate = df['flowrate'].tolist()
            else:
                dataframe = get_hydraulics(t)
                data_flowrate = dataframe['flowrate'].tolist()
        # return to OGS
        return (if_dyn_frate, data_flowrate)


# main
# initialize the tespy model of the bhe network
# load path of network model:
# loading the TESPy model
project_dir = os.getcwd()
print("Project dir is: ", project_dir)
nw = load_network('./pre/tespy_nw')
# set if print the network iteration info
nw.set_attr(iterinfo=False)

# create bhe dataframe of the network system from bhe_network.csv
df = create_dataframe()
n_BHE = np.size(df.iloc[:, 0])

# create local variables of the components label and connections label in
# network
localVars = locals()
data_index = df.index.tolist()
for i in range(n_BHE):
    for c in nw.conns.index:
        # bhe inlet and outlet conns
        if c.t.label == data_index[i]:  # inlet conns of bhe
            localVars['inlet_BHE' + str(i + 1)] = c
        if c.s.label == data_index[i]:  # outlet conns of bhe
            localVars['outlet_BHE' + str(i + 1)] = c

# time depended consumer thermal demand
if switch_dyn_demand == 'on':
    # import the name of bus from the network csv file
    bus_name = read_csv('./pre/tespy_nw/comps/bus.csv',
                    delimiter=';',
                    index_col=[0]).index[0]

# time depended flowrate
if switch_dyn_frate == 'on':
    # import the name of inlet connection from the network csv file
    inlet_name = read_csv('./pre/tespy_nw/conn.csv',
                    delimiter=';',
                    index_col=[0]).iloc[0,0]
    for c in nw.conns.index:
        # bhe inflow conns
        if c.s.label == inlet_name:  # inlet conns of bhe
            localVars['inlet_name'] = c

# instantiate BC objects referenced in OpenGeoSys
bc_bhe = BC()
