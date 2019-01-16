import sys
print(sys.version)
import os

import OpenGeoSys

from tespy import  cmp, con, nwk, hlp
from tespy import nwkr

import numpy as np
import pandas as pd
# %% create network dataframe
def create_dataframe():
    #return dataframe
    path_network =os.path.abspath('pre/bhe_network.csv')
    df_nw = pd.read_csv(path_network,delimiter=';', index_col=[0], dtype={'data_index':str})
    return(df_nw)
# %% TESPy hydraulic calculation process
def get_hydraulics():
    #refrigerant parameters
    refrig_density = 992.92 #kg/m3
    #solve imported network
    nw.solve(mode='design', design_file='pre/tespy_nw/results.csv', init_file='pre/tespy_nw/results.csv')
    #get flowrate #kg/s 
    for i in range(n_BHE):
        for c in nw.conns.index:
            if c.t.label == data_index[i]:#t:inlet comp, s:outlet comp
                df.loc[df.index[i],'flowrate'] = c.get_attr('m').val_SI
    #convert flowrate to velocity: #m^3/s
    for i in range(n_BHE):
        df.loc[df.index[i],'f_velocity'] = df.loc[df.index[i],'flowrate']/refrig_density
    return df
# %% TESPy Thermal calculation process
def get_thermal():
    #bhe network thermal parametrization
    for i in range(n_BHE):
        createVar['outlet_BHE'+ str(i+1)].set_attr( T= df.loc[data_index[i],'Tout_val'])

    # solving network
    nw.solve(mode='design')    

    #get Tin_val
    for i in range(n_BHE):
        df.loc[df.index[i],'Tin_val'] = createVar['inlet_BHE'+ str(i+1)].get_attr('T').val_SI

    return df['Tin_val'].tolist()
# %% OGS setting
# Dirichlet BCs
class BC(OpenGeoSys.BHENetwork):
    def initializeDataContainer(self):
        #convert dataframe to column list
        data_col_1 = df['Tin_val'].tolist()#'Tin_val'
        data_col_2 = df['Tout_val'].tolist()#'Tout_val'
        data_col_3 = df['Tout_node_id'].astype(int).tolist()#'Tout_node_id'

        return (True, data_col_1,data_col_2,data_col_3)
    def tespyThermalSolver(self, Tin_val, Tout_val):
        #read Tout_val to dataframe
        for i in range(n_BHE):
            df.loc[df.index[i],'Tout_val'] = Tout_val[i]
        #TESPy solver
        cur_cal_Tin_val = get_thermal()
        #check norm if network achieves the converge
        if_success = False
        pre_cal_Tin_val = Tin_val
        norm = np.linalg.norm(abs(np.asarray(pre_cal_Tin_val)-np.asarray(cur_cal_Tin_val)))
        if norm < 10e-6:
            if_success = True
        #return to OGS
        return (True, if_success, cur_cal_Tin_val)
    def tespyHydroSolver(self):
        dataframe = get_hydraulics()
        data_f_velocity = dataframe['f_velocity'].tolist()#'f_velocity'
        #return to OGS
        return (True, data_f_velocity)
# %% main
#initialize the tespy model of the bhe network
#load path of network model:
nw_path = os.path.abspath('pre/tespy_nw')
#load the model
nw = nwkr.load_nwk(nw_path)
#set if print the information of the network
nw.set_printoptions(print_level='none')

#create bhe dataframe of the network system from bhe_network.csv
df = create_dataframe()
n_BHE = np.size(df.iloc[:,0])
#bhes name
data_index = df.index.tolist()

#create global variables of the components label and connections label in network
createVar = locals() 

data_index = df.index.tolist()
for i in range(n_BHE):
    for c in nw.conns.index:
        #bhe inlet and outlet conns
        if c.t.label == data_index[i]:# inlet conns of bhe
            createVar['inlet_BHE'+ str(i+1)] = c
        if c.s.label == data_index[i]:# outlet conns of bhe
            createVar['outlet_BHE'+ str(i+1)] = c

# instantiate BC objects referenced in OpenGeoSys
bc_bhe = BC()