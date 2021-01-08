import mtest
import numpy as np
import matplotlib.pyplot as plt

m = mtest.MTest()
mtest.setVerboseMode(mtest.VerboseLevel.VERBOSE_QUIET)
m.setMaximumNumberOfSubSteps(20)
m.setModellingHypothesis("Axisymmetrical")
m.setBehaviour("generic", "src/libBehaviour.so", "ModCamClay_semiExplParaInit")

# Material constants (according to Modified Cam clay model Report)
nu = 0.3  # Poisson ratio
E = 2 * (1 + nu) * 20.0e6  # Young's modulus in Pa
la = 7.7e-2  # slope of the virgin consolidation line
ka = 6.6e-3  # slope of the swelling line
M = 1.2  # slope of the critical state line (CSL)
v0 = 1.788  # initial volume ratio
phi0 = 1 - 1 / v0  # Initial porosity
pc0 = 200.0e3  # Initial pre-consolidation pressure in Pa
pamb = 1.0e3  # Ambient pressure in Pa

# Loading programm
tMax = 1.0  # s , total time
nTime = 200
ltime = np.linspace(0.0, tMax, nTime)

p_axi = 587387  # axial pressure, +12614 for reaching CSL
p_con = 200000  # confining pressure
e_con = p_con * (1 - 2 * nu) / E
m.setImposedStress("SRR", {0: 0, 0.02: -p_con, 1.0: -p_con})
m.setImposedStress("STT", {0: 0, 0.02: -p_con, 1.0: -p_con})
# stress-controlled: works only until reaching the CSL
m.setImposedStress("SZZ", {0: 0, 0.02: -p_con, 1.0: -p_axi})
# strain-controlled: works, CSL reached asymptotically for EYY->inf
# m.setImposedStrain('EZZ', {0:0, 0.02:-e_con, 1.0:-130*e_con})

# Environment parameters
m.setExternalStateVariable("Temperature", 293.15)
m.setParameter("AmbientPressure", pamb)

# Material parameters
m.setMaterialProperty("YoungModulus", E)
m.setMaterialProperty("PoissonRatio", nu)
m.setMaterialProperty("CriticalStateLineSlope", M)
m.setMaterialProperty("SwellingLineSlope", ka)
m.setMaterialProperty("VirginConsolidationLineSlope", la)

# Initial values (only for the ParaInit code version!)
m.setMaterialProperty("InitialPreConsolidationPressure", pc0)
m.setMaterialProperty("InitialPorosity", phi0)

s = mtest.MTestCurrentState()
wk = mtest.MTestWorkSpace()
m.completeInitialisation()
m.initializeCurrentState(s)
m.initializeWorkSpace(wk)

# initialize output lists
pCurve = np.array([pamb])
qCurve = np.array([0.0])
eVCurve = np.array([0.0])
eQCurve = np.array([0.0])
lpCurve = np.array([0.0])
pcCurve = np.array([pc0])
phiCurve = np.array([phi0])
strains = np.empty(shape=(6, nTime))
stresses = np.empty(shape=(6, nTime))

# initialize yield functions
nPoints = 1000
pRange = np.empty(shape=(nTime, nPoints))
qFunct = np.empty(shape=(nTime, nPoints))

# run sim
for i in range(nTime - 1):
    m.execute(s, wk, ltime[i], ltime[i + 1])

    # output variables:
    pressure = pamb - (s.s1[0] + s.s1[1] + s.s1[2]) / 3
    vMstress = np.sqrt(
        0.5
        * (
            (s.s1[0] - s.s1[1]) ** 2
            + (s.s1[1] - s.s1[2]) ** 2
            + (s.s1[2] - s.s1[0]) ** 2
            + 3 * s.s1[3] ** 2
        )
    )
    epsilonV = s.e1[0] + s.e1[1] + s.e1[2]
    argument = (
        2
        * (
            s.e1[0] ** 2
            + s.e1[1] ** 2
            + s.e1[2] ** 2
            - epsilonV ** 2 / 3
            + 2 * s.e1[3] ** 2
        )
        / 3
    )
    vMstrain = np.sqrt(max(argument, 0))
    eplEquiv = s.getInternalStateVariableValue("EquivalentPlasticStrain")
    porosity = s.getInternalStateVariableValue("Porosity")
    pc = s.getInternalStateVariableValue("PreConsolidationPressure")
    eplV = s.getInternalStateVariableValue("PlasticVolumetricStrain")

    pCurve = np.append(pCurve, pressure)
    qCurve = np.append(qCurve, vMstress)
    phiCurve = np.append(phiCurve, porosity)
    pcCurve = np.append(pcCurve, pc)
    eVCurve = np.append(eVCurve, epsilonV)
    eQCurve = np.append(eQCurve, vMstrain)
    lpCurve = np.append(lpCurve, eplEquiv)

    for k in range(4):
        strains[k][i + 1] = s.e1[k]

    for k in range(4):
        stresses[k][i + 1] = s.s1[k]

    # calculate the yield surfaces in the p-q-space
    pRange[i] = np.linspace(0, pc, nPoints)
    qFunct[i] = M * np.sqrt(pRange[i] * (pc - pRange[i]))

# calculate the theoretical CSL
pRangeCSL = np.linspace(0, pc0 / 2, nPoints)
pRangeCSL = pRange[nTime - 2]
qFunctCSL = M * pRangeCSL

# calculate the analytical reference solution according to Peric (2006)
a = 3 * (1 - 2 * nu) / (2 * (1 + nu))
k = vMstress / (pressure - pc0)
c = (la - ka) / M
nP = 30
qRangeAna = np.linspace(0.0, vMstress, nP)
pRangeAna = np.linspace(pc0, pressure, nP)
x = qRangeAna / (M * pRangeAna)
y = 1 - qRangeAna / (k * pRangeAna)  # = (pc0/pRangeAna)

v0xEpsQp = np.log(
    np.power((1 - x), c * k / (M - k)) * np.power((1 + x), c * k / (M + k))
) - 2 * c * np.arctan(x)
v0xEpsQe = np.log(np.power(y, 2 * c / (k / M - M / k) - ka * k / (3 * a)))
v0xEpsQ = v0xEpsQe + v0xEpsQp

# print some final values
print("Triaxial test:")
print("final normal strain in z direction: ", s.e1[1])
print("final normal stress in z direction: ", s.s1[1], "Pa")
print("final von Mises stress: ", vMstress, "Pa")
print("final hydrostatic pressure: ", pressure, "Pa")
print("final pre-consolidation pressure: ", pc, "Pa")
print("confining strain in z direction: ", e_con)

# plots
fig, ax = plt.subplots()
ax.set_title("Numerical solution versus analytical solution (Peric, 2006)")
ax.scatter(v0xEpsQ / v0, qRangeAna / 1e3, label="analytical")
ax.plot(eQCurve, qCurve / 1e3, color="black", label="numerical")
ax.set_xlabel("$\epsilon_{q}$")
ax.set_ylabel("q / kPa")
ax.grid()
ax.legend()
fig.savefig("out/ModCamClay_TriaxStudy_NumVsAnal.pdf")

fig, ax = plt.subplots()
ax.set_title("Loading trajectories in the stress space")
ax.plot(pCurve, qCurve, label="$(p,q)$ / Pa")
plt.quiver(
    pCurve[:-1],
    qCurve[:-1],
    pCurve[1:] - pCurve[:-1],
    qCurve[1:] - qCurve[:-1],
    scale_units="xy",
    angles="xy",
    scale=1,
)
ax.plot(pRangeCSL, qFunctCSL, label="CSL")
for k in range(0, nTime - 1, 20):
    ax.plot(pRange[k], qFunct[k])
ax.plot(pRange[nTime - 2], qFunct[nTime - 2])
ax.set_xlabel("$p$ / Pa")
ax.set_ylabel("$q$ / Pa")
ax.grid()
ax.legend()
fig.savefig("out/ModCamClay_TriaxStudy_YieldSurface.pdf")

fig, ax = plt.subplots()
ax.plot(ltime, pCurve, label="$p$ / Pa")
ax.plot(ltime, qCurve, label="$q$ / Pa")
ax.plot(ltime, pcCurve, label="$pc$ / Pa")
ax.set_xlabel("$t$ / s")
ax.set_ylabel("stress / Pa")
ax.grid()
ax.legend()

fig, ax = plt.subplots()
ax.plot(ltime, stresses[0][:], label="$\sigma_{rr}$ / MPa")
ax.plot(ltime, stresses[2][:], label="$\sigma_{\phi\phi}$ / MPa")
ax.plot(ltime, stresses[1][:], label="$\sigma_{zz}$ / MPa")
ax.plot(ltime, stresses[3][:], label="$\sigma_{rz}$ / MPa")
ax.set_xlabel("$t$ / s")
ax.set_ylabel("stress / MPa")
ax.grid()
ax.legend()

fig, ax = plt.subplots()
ax.plot(ltime, strains[0][:], color="red", label="$\epsilon_{rr}$")
ax.plot(ltime, strains[2][:], "--", color="green", label="$\epsilon_{\phi\phi}$")
ax.plot(ltime, strains[1][:], label="$\epsilon_{zz}$")
ax.plot(ltime, strains[3][:], color="black", label="$\epsilon_{rz}$")
ax.set_xlabel("$t$ / s")
ax.set_ylabel("strain")
ax.grid()
ax.legend(loc="lower left")
fig.savefig("out/ModCamClay_TriaxStudy_Strains.pdf")

fig, ax = plt.subplots()
ax.plot(ltime, phiCurve - phi0, label="$\phi-\phi_0$")
ax.plot(ltime, eVCurve, label="$\epsilon_{V}$")
ax.plot(ltime, lpCurve, label="$\epsilon_{eq}$")
ax.set_xlabel("$t$ / s")
ax.set_ylabel(" ")
ax.grid()
ax.legend()

plt.show()
