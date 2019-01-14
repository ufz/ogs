# -*- coding: utf-8 -*-
"""
Created on Wed Jun 13 11:21:54 2018

@author: geche
object:
    preprocess for prepare a python script for buliding a muti bhe network draft
"""
#################################
import numpy as np
import scipy
import os
import math
import pandas as pd

##############################################################################
#file SETTINGS:
#########################
#os.mkdir(""C://george//PhD//UFZ//Prof. Shao//Übung//31//tespy_impl") 
os.chdir("C://george//PhD//UFZ//Prof. Shao//Übung//31//TESPY_sc//tespy_new_branch_sc//testszenario//case4_lagrebhearray_subbtes") 
output_file_path = os.getcwd() 
#input_file_path = "C://george//PhD//UFZ//Prof. Shao//Übung//31//tespy_impl//???"

#%% user default setting
# splitter and merge parameters
sub_sp_num = 5

#%% user default setting
# pipeline parameters
pipe_typ =1 
#pipe parameters: length [m], diameter [m], roughness[-]
pipe = pd.DataFrame(columns=['l', 'd', 'r'])
pipe.loc[0] =[100,0.1,100]
pipe.loc[1] = [100,0.2,300]

bhe = pd.DataFrame(columns=['l', 'd', 'r'])
bhe.loc[0] =[50,0.02733,100]

#how many bhes in one bhe array?
bhes_num = 5

#%%text plot
#components:
with open('bhe_array.tec', 'w') as f:
##spitter and merge
    #components:
    for i in range(1,sub_sp_num+1):
        str_input = f'\nsp{i} = cmp.splitter(\'splitter{i}\',num_out = 5)'
        f.write(str_input)
    for i in range(1,sub_sp_num+1):
        str_input = f'\nmg{i} = cmp.merge(\'merge{i}\',num_in = 5)'
        f.write(str_input)
    #connections:
    for i in range(1,sub_sp_num+1):
        str_input = f'\nsp_sp{i} = con.connection(sp, \'out{i}\', sp{i}, \'in1\')'
        f.write(str_input)
    for i in range(1,sub_sp_num+1):
        str_input = f'\nmg{i}_mg = con.connection(mg{i}, \'out1\', mg, \'in{i}\')'
        f.write(str_input)
    #splitter,bhes,merge:
    for i in range(1,sub_sp_num+1):
        str_input = f'\nsp{i}_bhe{i} = con.connection(sp{i}, \'out1\', bhe{i}.inlet, \'in1\')'
        f.write(str_input)
    for i in range(1,sub_sp_num+1):
        str_input = f'\nbhe{i}_mg{i} = con.connection(bhe{i}.outlet, \'out1\', mg{i}, \'in1\')'
        f.write(str_input)
    #add connections:
    for i in range(1,sub_sp_num+1):
        str_input = f'\nsp_sp{i},sp{i}_bhe{i},\n\
bhe{i}_mg{i},mg{i}_mg,'
        f.write(str_input)
    f.write('\n')
    
##pipelines:
    for i in range(bhes_num):
        str_input = f'\nbhe{i+1}.set_attr(dT_pb{i}=0.0,L_pb{i}={pipe.iloc[0,0]},ks_pb{i}={pipe.iloc[0,2]},D_pb{i}={pipe.iloc[0,1]},\n\
dT_pf{i}=0.0,L_pf{i}={pipe.iloc[0,0]},ks_pf{i}={pipe.iloc[0,2]},D_pf{i}={pipe.iloc[0,1]},\n\
Q{i}=1000,L_bhe{i}={bhe.iloc[0,0]},ks_bhe{i}={bhe.iloc[0,2]},D_bhe{i}={bhe.iloc[0,1]},\n\
hydro_group =\'HW\')'
        f.write(str_input)
    for i in range(bhes_num):
        str_input = f'\nbhe{i+1}.set_attr(dT_pb{i}=0.0,L_pb{i}={pipe.iloc[1,0]},ks_pb{i}={pipe.iloc[1,2]},D_pb{i}={pipe.iloc[1,1]},\n\
dT_pf{i}=0.0,L_pf{i}={pipe.iloc[1,0]},ks_pf{i}={pipe.iloc[1,2]},D_pf{i}={pipe.iloc[1,1]},\n\
Q{i}=1000,L_bhe{i}={bhe.iloc[0,0]},ks_bhe{i}={bhe.iloc[0,2]},D_bhe{i}={bhe.iloc[0,1]},\n\
hydro_group =\'HW\')'
        f.write(str_input)
f.close
