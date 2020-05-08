#!/ usr / bin / env python
#coding : utf - 8

#![Logo TUBAF]( \
    https: //tu-freiberg.de/sites/default/files/media/freiberger-alumni-netzwerk-6127/wbm_orig_rgb_0.jpg)
#
#Run skript for Mohr Coulomb with anisotropy via stress scaling or \
    embedded weakness planes.
#Comments to:
#
#* Prof.Dr.Thomas Nagel
#Chair of Soil Mechanics and Foundation Engineering
#Geotechnical Institute
#Technische Universität Bergakademie Freiberg.*
#
#https:  // tu-freiberg.de/en/fakultaet3/gt/soilmechanics
#
#Tests and parameter sets based on:
#
#Malcom, M.A.M.(2018).Analysis of underground excavations in argillaceous hard \
                 soils -                                                       \
             weak rocks.Technical University of Catalonia.
#
#Ismael, M., &Konietzky,                                                   \
    H.(2017)                                                               \
        .Integration of Elastic Stiffness Anisotropy into Ubiquitous Joint \
            Model.Procedia Engineering,                                    \
    191, 1032–1039. https:  // doi.org/10.1016/j.proeng.2017.05.276

#In[1]:

import numpy as np import os as os import pyvista as pv import matplotlib.pyplot as plt

#In[2]:

#Some plot settings
    plt.style.use('seaborn-deep') plt.rcParams['lines.linewidth'] = 2.0 plt.rcParams['lines.color'] = 'black' plt.rcParams['legend.frameon'] = True plt.rcParams['font.family'] = 'serif' plt.rcParams['legend.fontsize'] = 14 plt.rcParams['font.size'] = 14 plt.rcParams['axes.spines.right'] = False plt.rcParams['axes.spines.top'] = False plt.rcParams['axes.spines.left'] = True plt.rcParams['axes.spines.bottom'] = True plt.rcParams['axes.axisbelow'] = True plt.rcParams['figure.figsize'] =(12, 6)

#In[3]:

                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 cmd = 'rm triax*.pvd triax_*pcs_0*.vtu triax_*1e0_*.prj' os.system(cmd)

#In[4]:

                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        beta = np.linspace(0, np.pi / 2, 28)

#In[5]:

                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               for i in beta : print("Generating input for %i °" %(int(np.round(i * 180 / np.pi, 0)))) filename = "triax_1e0_" + str(int(np.round(i * 180. / np.pi, 0))) + ".prj" cmd = "sed 's/<values>1 0/<values>" + str(np.cos(i)) + " " + str(np.sin(i)) + "/g' triax_original.prj > " + filename os.system(cmd) cmd = "sed 's/<values>0 1/<values>" + str(- np.sin(i)) + " " + str(np.cos(i)) + "/g' " + filename + " -i" os.system(cmd) cmd = "sed 's/<prefix>triax_0/<prefix>triax_" + str(int(np.round(i * 180. / np.pi, 0))) + "/g' " + filename + " -i" os.system(cmd) print("Running simulation for %i °" %(int(np.round(i * 180 / np.pi, 0)))) cmd = "~/ogs_release/bin/ogs " + filename os.system(cmd)

#In[6]:

                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    UCS = beta * 0. for n, i in enumerate(beta) :
#star, as ts value unknown due to adaptive time stepping
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                a = 'triax_' + str(int(np.round(i * 180. / np.pi, 0))) + '_pcs_0_ts_*_t_2.000000.vtu' filename = get_ipython().getoutput('ls $a') data = pv.read(filename[0]) sigma = data.get_array('sigma') UCS[n] = - sigma[0][1]

#In[7]:

                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  fig, ax = plt.subplots() ax.plot(90 - beta * 180 / np.pi, UCS / 1e3, marker = 'd') ax.set_xlabel('Angle of fracture plane / °') ax.set_ylabel('UCS / kPa') ax.set_ylim(2, 10) ax.set_xlim(0, 90);
ax.set_title('isotropic elasticity')
fig.savefig('orientation_dependent_strength.pdf')

#In[8]:


for i in beta:
    print("Generating input for %i °" %(int(np.round(i*180/np.pi,0))))
    filename = "triax_ortho_1e0_"+str(int(np.round(i*180./np.pi,0)))+".prj"
    cmd = "sed 's/<values>1 0/<values>" + str(np.cos(i)) + " " + str(np.sin(i))+ "/g' triax_ORTHO_original.prj > " + filename
    os.system(cmd)
    cmd = "sed 's/<values>0 1/<values>" + str(-np.sin(i)) + " " + str(np.cos(i))+ "/g' " + filename + " -i"
    os.system(cmd)
    cmd = "sed 's/<prefix>triax_ortho_0/<prefix>triax_ortho_" +str(int(np.round(i*180./np.pi,0))) + "/g' " + filename + " -i"
    os.system(cmd)
    print("Running simulation for %i °" %(int(np.round(i*180/np.pi,0))))
    cmd = "~/ogs_release/bin/ogs " + filename
    os.system(cmd)

#In[9]:


UCS = beta * 0.
for n,i in enumerate(beta):
#star, as ts value unknown due to adaptive time stepping
    a = 'triax_ortho_'+str(int(np.round(i*180./np.pi,0)))+'_pcs_0_ts_*_t_2.000000.vtu'
    filename = get_ipython().getoutput('ls $a')
    data = pv.read(filename[0])
    sigma = data.get_array('sigma')
    UCS[n] = -sigma[0][1]

#In[10]:


fig, ax = plt.subplots()
ax.plot(90 - beta*180/np.pi,UCS/1e6,marker='d')
ax.set_xlabel('Angle of fracture plane / °')
ax.set_ylabel('UCS / MPa')
ax.set_ylim(0,120)
ax.set_xlim(0,90);
ax.set_title('transversely isotropic elasticity')
fig.savefig('orientation_dependent_strength_ortho.pdf')

#In[11]:


for i in beta:
    print("Generating input for %i °" %(int(np.round(i*180/np.pi,0))))
    filename = "triax_aniso_1e0_"+str(int(np.round(i*180./np.pi,0)))+".prj"
    cmd = "sed 's/<values>1 0/<values>" + str(np.cos(i)) + " " + str(np.sin(i))+ "/g' triax_aniso_original.prj > " + filename
    os.system(cmd)
    cmd = "sed 's/<values>0 1/<values>" + str(-np.sin(i)) + " " + str(np.cos(i))+ "/g' " + filename + " -i"
    os.system(cmd)
    cmd = "sed 's/<prefix>triax_aniso_0/<prefix>triax_aniso_" +str(int(np.round(i*180./np.pi,0))) + "/g' " + filename + " -i"
    os.system(cmd)
    print("Running simulation for %i °" %(int(np.round(i*180/np.pi,0))))
    cmd = "~/ogs_release/bin/ogs " + filename
    os.system(cmd)

#In[12]:


UCS = beta * 0.
for n,i in enumerate(beta):
#star, as ts value unknown due to adaptive time stepping
    a = 'triax_aniso_'+str(int(np.round(i*180./np.pi,0)))+'_pcs_0_ts_*_t_2.000000.vtu'
    filename = get_ipython().getoutput('ls $a')
    data = pv.read(filename[0])
    sigma = data.get_array('sigma')
    UCS[n] = -sigma[0][1]

#In[13]:


fig, ax = plt.subplots()
ax.plot(90 - beta*180/np.pi,UCS/1e6,marker='d')
ax.set_xlabel('Angle of fracture plane / °')
ax.set_ylabel('UCS / MPa')
ax.set_ylim(0,120)
ax.set_xlim(0,90);
ax.set_title('transversely isotropic elasticity') fig.savefig(
    'orientation_dependent_strength_aniso.pdf')

#In[]:
