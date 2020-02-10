###
# Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
# Distributed under a Modified BSD License.
# See accompanying file LICENSE.txt or
# http://www.opengeosys.org/project/license
###

# Execute this file to generate TESPy network csv files
from tespy.networks import network
from tespy.components import (sink, source, splitter, merge,
                              pump, heat_exchanger_simple)
from tespy.connections import connection, ref, bus
from tespy.tools.characteristics import char_line
from tespy.tools.data_containers import dc_cc

import numpy as np

# %% network
btes = network(fluids=['water'],
                   T_unit='K',
                   p_unit='bar',
                   h_unit='kJ / kg',
                   T_range=[273.25, 373.15],
                   p_range=[1, 20],
                   h_range=[1, 1000])

# %% components
fc_in = source('from consumer inflow')
fc_out = sink('from consumer outflow')

pu = pump('pump')

sp = splitter('splitter', num_out=3)

# bhe:
bhe_name = 'BHE1'
assert 'BHE1' in bhe_name, "BHE should be named with 'BHE1'"
bhe1 = heat_exchanger_simple(bhe_name)
bhe_name = 'BHE2'
assert 'BHE2' in bhe_name, "BHE should be named with 'BHE2'"
bhe2 = heat_exchanger_simple(bhe_name)
bhe_name = 'BHE3'
assert 'BHE3' in bhe_name, "BHE should be named with 'BHE3'"
bhe3 = heat_exchanger_simple(bhe_name)

mg = merge('merge', num_in=3)

cons = heat_exchanger_simple('consumer')

# %% connections
fc_pu = connection(fc_in, 'out1', pu, 'in1')

pu_sp = connection(pu, 'out1', sp, 'in1')

sp_bhe1 = connection(sp, 'out1', bhe1, 'in1')
sp_bhe2 = connection(sp, 'out2', bhe2, 'in1')
sp_bhe3 = connection(sp, 'out3', bhe3, 'in1')

bhe1_mg = connection(bhe1, 'out1', mg, 'in1')
bhe2_mg = connection(bhe2, 'out1', mg, 'in2')
bhe3_mg = connection(bhe3, 'out1', mg, 'in3')

mg_cons = connection(mg, 'out1', cons, 'in1')

cons_fc = connection(cons, 'out1', fc_out, 'in1')

btes.add_conns(fc_pu, pu_sp, sp_bhe1, sp_bhe2, sp_bhe3, bhe1_mg, bhe2_mg,
               bhe3_mg, mg_cons, cons_fc)


# %% paramerization
## components paramerization
# pump
# flow_char
# provide volumetric flow in m^3 / s
x = np.array([
    0.00, 0.00001952885971862, 0.00390577194372, 0.005858657915586,
    0.007811543887448, 0.00976442985931, 0.011717315831173, 0.013670201803035,
    0.015623087774897, 0.017575973746759, 0.019528859718621, 0.021481745690483,
    0.023434631662345, 0.025387517634207, 0.027340403606069, 0.029293289577931,
    0.031246175549793, 0.033199061521655, 0.035151947493517, 0.037104833465379,
    0.039057719437241, 0.041010605409104, 0.042963491380966, 0.044916377352828,
    0.04686926332469, 0.048822149296552, 0.050775035268414, 0.052727921240276,
    0.054680807212138, 0.056633693184
])

# provide head in Pa
y = np.array([
    0.47782539, 0.47725723, 0.47555274, 0.47271192, 0.46873478, 0.46362130,
    0.45737151, 0.44998538, 0.44146293, 0.43180416, 0.4220905, 0.40907762,
    0.39600986, 0.38180578, 0.36646537, 0.34998863, 0.33237557, 0.31362618,
    0.29374046, 0.27271841, 0.25056004, 0.22726535, 0.20283432, 0.17726697,
    0.15056329, 0.12272329, 0.09374696, 0.06363430, 0.03238531, 0.00000000
]) * 1e5

char = char_line(x=x, y=y)
pu.set_attr(flow_char=dc_cc(func=char, is_set=True))
pu.set_attr(eta_s=0.90)

# bhes
bhe1.set_attr(D=0.013665, L=100, ks=0.00001)
bhe2.set_attr(D=0.013665, L=100, ks=0.00001)
bhe3.set_attr(D=0.013665, L=100, ks=0.00001)

# consumer
cons.set_attr(D=0.2, L=20, ks=0.00001)
# busses
heat = bus('consumer heat demand')
heat.add_comps({'c': cons, 'p': 'P'})
btes.add_busses(heat)
# consumer heat demand
heat.set_attr(P=-3000)  # W


## connection parametrization
# system inlet
inflow_head = 2  # bar

fc_pu.set_attr(p=inflow_head, m=0.6, fluid={'water': 1})

# for BHEs:
# Tout:
bhe1_mg.set_attr(T=303.15)
bhe2_mg.set_attr(T=303.15)
bhe3_mg.set_attr(T=303.15)

# imposed boundary condition: ensure all heat from BHEs are consumed on 'consumer'
pu_sp.set_attr(h=ref(cons_fc, 1, 0))

# %% solve
btes.solve('design')
#btes.print_results()

# %% save to csv:
btes.save('tespy_nw', structure=True)
