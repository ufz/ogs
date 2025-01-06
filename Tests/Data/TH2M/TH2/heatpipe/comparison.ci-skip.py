# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.14.5
#   kernelspec:
#     display_name: Python 3.8.10 64-bit
#     metadata:
#       interpreter:
#         hash: 5b3ded1ccb95c1d9bd405e7b823d9e85424cde40fbb5985eb47e999ef50e15b4
#     name: python3
# ---

# %%
import matplotlib.pyplot as plt
import pandas as pd
import vtuIO
import numpy as np

# %%
analytical_solution = pd.read_csv("analytical.csv")
numerical_solution = vtuIO.PVDIO("results_heatpipe_rough.pvd", dim=2)

# %%
time = 1e7
x = np.linspace(0.0, 1.0, 201)
line = [(i, 0.0025, 0) for i in x]

T_num = numerical_solution.read_point_set_data(time, "temperature", pointsetarray=line)
pGR_num = numerical_solution.read_point_set_data(
    time, "gas_pressure", pointsetarray=line
)
pCap_num = numerical_solution.read_point_set_data(
    time, "capillary_pressure", pointsetarray=line
)
xnCG_num = numerical_solution.read_point_set_data(time, "xnCG", pointsetarray=line)
sL_num = numerical_solution.read_point_set_data(time, "saturation", pointsetarray=line)

# %%
analytical_solution

# %%
z_a = analytical_solution["z"].to_numpy()
T_a = analytical_solution["T"].to_numpy()
pGR_a = analytical_solution["pGR"].to_numpy()
pCap_a = analytical_solution["pCap"].to_numpy()
xnCG_a = analytical_solution["xnCG"].to_numpy()
sL_a = analytical_solution["sL"].to_numpy()

# %%
fig, ax = plt.subplots(figsize=(12, 6))

ax.plot(x, T_num, "*", label="T_numerical")
ax.plot(z_a, T_a, "-", label="T_analytical")

plt.xlim([0, 1])
plt.legend(loc="upper right", bbox_to_anchor=(0.9, 0.3))
plt.show()

# %%
fig, ax = plt.subplots(figsize=(12, 6))

ax.plot(x, sL_num, "*", label="sL_numerical")
ax.plot(z_a, sL_a, "-", label="sL_analytical")

ax.plot(x, xnCG_num, "*", label="xnCG_numerical")
ax.plot(z_a, xnCG_a, "-", label="xnCG_analytical")

plt.xlim([0, 1])
plt.legend(loc="upper right", bbox_to_anchor=(0.9, 0.8))
plt.show()

# %%
fig, ax = plt.subplots(figsize=(12, 6))

ax.plot(x, pGR_num, "*", label="pGR_numerical")
ax.plot(z_a, pGR_a, "-", label="pGR_analytical")

ax.plot(x, pCap_num, "*", label="pCap_numerical")
ax.plot(z_a, pCap_a, "-", label="pCap_analytical")

plt.xlim([0, 1])
plt.ylim([0, 180000])

plt.legend(loc="upper right", bbox_to_anchor=(0.9, 0.8))
plt.show()

# %%
