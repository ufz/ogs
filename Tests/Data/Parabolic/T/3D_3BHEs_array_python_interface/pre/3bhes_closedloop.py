# SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
# SPDX-License-Identifier: BSD-3-Clause

# Execute this file to generate TESPy network JSON file
import numpy as np
from tespy.components import (
    CycleCloser,
    Merge,
    Pump,
    SimpleHeatExchanger,
    Splitter,
)
from tespy.connections import Connection
from tespy.networks import Network
from tespy.tools.characteristics import CharLine
from tespy.tools.fluid_properties.wrappers import IAPWSWrapper

# %% network
btes = Network()
btes.units.set_defaults(
    temperature="K", pressure="bar", pressure_difference="bar", enthalpy="kJ / kg"
)

# %% components
fc = CycleCloser("cycle closer")
pu = Pump("pump")
sp = Splitter("splitter", num_out=3)

# bhe:
bhe1 = SimpleHeatExchanger("BHE1")
bhe2 = SimpleHeatExchanger("BHE2")
bhe3 = SimpleHeatExchanger("BHE3")

mg = Merge("merge", num_in=3)
cons = SimpleHeatExchanger("consumer")

## components paramerization
# pump
# flow_char
# provide volumetric flow in m^3 / s
x = np.array(
    [
        0.00,
        0.00001952885971862,
        0.00390577194372,
        0.005858657915586,
        0.007811543887448,
        0.00976442985931,
        0.011717315831173,
        0.013670201803035,
        0.015623087774897,
        0.017575973746759,
        0.019528859718621,
        0.021481745690483,
        0.023434631662345,
        0.025387517634207,
        0.027340403606069,
        0.029293289577931,
        0.031246175549793,
        0.033199061521655,
        0.035151947493517,
        0.037104833465379,
        0.039057719437241,
        0.041010605409104,
        0.042963491380966,
        0.044916377352828,
        0.04686926332469,
        0.048822149296552,
        0.050775035268414,
        0.052727921240276,
        0.054680807212138,
        0.056633693184,
    ]
)

# provide head in Pa
y = (
    np.array(
        [
            0.47782539,
            0.47725723,
            0.47555274,
            0.47271192,
            0.46873478,
            0.46362130,
            0.45737151,
            0.44998538,
            0.44146293,
            0.43180416,
            0.4220905,
            0.40907762,
            0.39600986,
            0.38180578,
            0.36646537,
            0.34998863,
            0.33237557,
            0.31362618,
            0.29374046,
            0.27271841,
            0.25056004,
            0.22726535,
            0.20283432,
            0.17726697,
            0.15056329,
            0.12272329,
            0.09374696,
            0.06363430,
            0.03238531,
            0.00000000,
        ]
    )
    * 1e5
)
char = CharLine(x=x, y=y)
pu.set_attr(flow_char={"char_func": char, "is_set": True})
pu.set_attr(eta_s=0.90)

# bhes
bhe1.set_attr(D=0.013665, L=100, ks=0.00001)
bhe2.set_attr(D=0.013665, L=100, ks=0.00001)
bhe3.set_attr(D=0.013665, L=100, ks=0.00001)

# consumer
cons.set_attr(pr=1)
# consumer heat demand
cons.set_attr(Q=-3000)  # W

# %% connections
fc_pu = Connection(fc, "out1", pu, "in1", label="fc_pu")
pu_sp = Connection(pu, "out1", sp, "in1", label="pu_sp")

sp_bhe1 = Connection(sp, "out1", bhe1, "in1", label="sp_bhe1")
sp_bhe2 = Connection(sp, "out2", bhe2, "in1", label="sp_bhe2")
sp_bhe3 = Connection(sp, "out3", bhe3, "in1", label="sp_bhe3")

bhe1_mg = Connection(bhe1, "out1", mg, "in1", label="bhe1_mg")
bhe2_mg = Connection(bhe2, "out1", mg, "in2", label="bhe2_mg")
bhe3_mg = Connection(bhe3, "out1", mg, "in3", label="bhe3_mg")

mg_cons = Connection(mg, "out1", cons, "in1", label="mg_cons")
cons_fc = Connection(cons, "out1", fc, "in1", label="cons_fc")

btes.add_conns(
    fc_pu, pu_sp, sp_bhe1, sp_bhe2, sp_bhe3, bhe1_mg, bhe2_mg, bhe3_mg, mg_cons, cons_fc
)

## connection parametrization
# system inlet
# The IAPWS-IF97 engine is used instead of TESPy's CoolProp default: CoolProp
# crashes inside OGS's embedded Python, see bhe_tespy_compat.py. The engine is
# propagated from here to all other connections of the network.
fc_pu.set_attr(p=2, fluid={"water": 1}, fluid_engines={"water": IAPWSWrapper})

# for BHEs:
# Tout:
bhe1_mg.set_attr(T=303.15)
bhe2_mg.set_attr(T=303.15)
bhe3_mg.set_attr(T=303.15)


# %% solve
btes.solve("design")
# btes.print_results()

# %% save to json:
btes.export("tespy_nw_closedloop.json")
