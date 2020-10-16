import mtest
import numpy as np
import matplotlib.pyplot as plt


GM0 = 9.54e3
KM0 = 2.78e4
GK0 = 6.27e4
mK = -0.254
mvK = -0.327
mvM = -0.267
etaK0 = 1.66e5
etaM0 = 4.03e7
tmax = 15.0


m = mtest.MTest()
mtest.setVerboseMode(mtest.VerboseLevel.VERBOSE_QUIET)
m.setMaximumNumberOfSubSteps(10)
m.setBehaviour("generic", "src/libBehaviour.so", "Lubby2")
m.setExternalStateVariable("Temperature", 293.15)
m.setImposedStress("SXX", 0.0)
m.setImposedStress("SYY", 0.0)
m.setImposedStress("SZZ", 0.0)
m.setImposedStress("SXZ", 0.0)
m.setImposedStress("SYZ", 0.0)
m.setImposedStress(
    "SXY", {0: 0, 1.0e-5: 5 * np.sqrt(2.0), 15: 5 * np.sqrt(2.0)}
)  # Kelvin Mapping of shear components!!

m.setMaterialProperty(
    "YoungModulus", (9 * GM0 * KM0) / (3 * KM0 + GM0)
)  # this was the wrong equation
m.setMaterialProperty(
    "PoissonRatio", (3 * KM0 - 2 * GM0) / (2 * GM0 + 6 * KM0)
)  # this was the wrong equation
m.setMaterialProperty("KelvinShearModulus", GK0)
m.setMaterialProperty("KelvinViscosity", etaK0)
m.setMaterialProperty("MaxwellViscosity", etaM0)
m.setMaterialProperty("KelvinElasticParameter", mK)
m.setMaterialProperty("MaxwellViscoParameter", mvM)  # these were switched
m.setMaterialProperty("KelvinViscoParameter", mvK)  # these were switched

sig_xy = 5.0
sig_eff = np.sqrt(3.0) * sig_xy

GK = GK0 * np.exp(mK * sig_eff)
etaK = etaK0 * np.exp(mvK * sig_eff)
etaM = etaM0 * np.exp(mvM * sig_eff)


eps_xy = (
    lambda t: (
        (1.0 / GM0 + t / etaM) * sig_xy
        + 1.0 / GK * (1.0 - np.exp(-GK / etaK * t)) * sig_xy
    )
    / 2.0
)

s = mtest.MTestCurrentState()
wk = mtest.MTestWorkSpace()
m.completeInitialisation()
m.initializeCurrentState(s)
m.initializeWorkSpace(wk)

ltime = np.append(np.linspace(0, 1e-5, 10), np.linspace(2.0e-5, tmax, 200))

numerical = np.array([0.0])

# run sim
for i in range(len(ltime) - 1):
    m.execute(s, wk, ltime[i], ltime[i + 1])
    numerical = np.append(
        numerical, s.e1[3] / np.sqrt(2.0)
    )  # Kelvin mapping backwards!

fig, ax = plt.subplots()
ax.plot(ltime, numerical * 100.0, label="numerical")
ax.plot(ltime, eps_xy(ltime) * 100.0, label="analytical")
ax.set_xlabel("$t$ / d")
ax.set_ylabel("$\\epsilon_{xy}$ / %")
fig.legend()
fig.savefig("lubby_mfront.pdf")
