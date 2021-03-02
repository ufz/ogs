import mtest
import numpy as np
import matplotlib.pyplot as plt


# Material constants
E = 150.0e3  # Young's modulus
nu = 0.3  # Poisson ratio
la = 0.0077  # slope of the virgin consolidation line
ka = 0.00066  # slope of the swelling line
M = 1.5  # slope of the critical state line (CSL)
v0 = 1.788  # initial volume ratio
pamb = 1e-3  # Ambient pressure
pc0 = 30  # Initial pre-consolidation pressure
phi0 = 1 - 1 / v0  # Initial porosity

# Loading constants
p = pc0 / 2  # Mpa, prescribed hydrostatic pressure
exy = 1e-2  # -  , prescribed shear strain
tMax = 1.0  # s  , total time

variableParameters = {
    1: "p",
    2: "M",
    3: "pc0",
}

# Choose parameter to be varied in this study
variedParameter = variableParameters[1]

# Lists for parameter effect studies
if variedParameter == "p":
    valueList = [pc0 / 4, pc0 / 2, 3 * pc0 / 4]  # typical range
    # valueList = [0,pc0/4,pc0/2,3*pc0/4,pc0,5*pc0/4]   # total range incl. 0
    # valueList = [0,pc0/10,pc0/6,pc0/3,pc0/2]          # softening range
    # valueList = [pc0,5*pc0/6,2*pc0/3,pc0/2]           # hardening range
    # valueList = [pc0,pc0+10,pc0+20,pc0+30]            # beyond pc
if variedParameter == "M":
    valueList = [M / 10, M / 2, M, 2 * M]
if variedParameter == "pc0":
    valueList = [pc0 / 4, 2 * pc0]

runs = max(len(valueList), 1)
ltime = np.linspace(0.0, tMax, 500)
nTime = len(ltime)
nOutput = 8
results = np.empty(shape=(nOutput, runs, nTime))
results[:][:][0] = 0.0

pEndvalues = np.empty(runs)
qEndvalues = np.empty(runs)
nPoints = 20
pRangeFinal = np.empty(shape=(runs, nPoints))
qFunctFinal = np.empty(shape=(runs, nPoints))

for k in range(runs):
    print("Materialset:", k)

    if variedParameter == "p":
        p = valueList[k]
        prelabel = "p ="
    if variedParameter == "M":
        M = valueList[k]
        prelabel = "$M=$"
    if variedParameter == "pc0":
        pc0 = valueList[k]
        prelabel = "$p_{c0}=$"

    m = mtest.MTest()
    mtest.setVerboseMode(mtest.VerboseLevel.VERBOSE_QUIET)
    m.setMaximumNumberOfSubSteps(10)
    m.setBehaviour("generic", "src/libBehaviour.so", "ModCamClay_semiExplParaInit")

    m.setExternalStateVariable("Temperature", 293.15)
    m.setImposedStress("SXX", {0: 0, 0.5: -p})
    m.setImposedStress("SYY", {0: 0, 0.5: -p})
    m.setImposedStress("SZZ", {0: 0, 0.5: -p})
    m.setImposedStrain("EXY", {0: 0, 0.5: 0, 1: exy})

    m.setParameter("AmbientPressure", pamb)

    m.setMaterialProperty("YoungModulus", E)
    m.setMaterialProperty("PoissonRatio", nu)
    m.setMaterialProperty("CriticalStateLineSlope", M)
    m.setMaterialProperty("SwellingLineSlope", ka)
    m.setMaterialProperty("VirginConsolidationLineSlope", la)
    # only for the ParaInit code version!
    m.setMaterialProperty("InitialPreConsolidationPressure", pc0)
    m.setMaterialProperty("InitialPorosity", phi0)

    m.setInternalStateVariableInitialValue("PreConsolidationPressure", pc0)
    m.setInternalStateVariableInitialValue("Porosity", phi0)
    m.setInternalStateVariableInitialValue("VolumeRatio", v0)

    s = mtest.MTestCurrentState()
    wk = mtest.MTestWorkSpace()
    m.completeInitialisation()
    m.initializeCurrentState(s)
    m.initializeWorkSpace(wk)

    # initialize output lists (other than zero)
    results[2][k][0] = v0
    results[3][k][0] = phi0
    results[6][k][0] = pc0

    # calculate the initial yield surface in the p-q-space
    pRange = np.linspace(0, pc0, 100)  # + pamb
    qFunct = M * np.sqrt((pc0 - pRange) * pRange)

    # run simulation
    for i in range(nTime - 1):
        m.execute(s, wk, ltime[i], ltime[i + 1])

        # output variables:
        # s0, e0: stress/strain at time[i]
        # s1, e1: stress/strain at time[i+1]
        pressure = -(s.s1[0] + s.s1[1] + s.s1[2]) / 3
        vMstress = np.sqrt(1.5) * s.s1[3]
        epsilonV = s.e1[0] + s.e1[1] + s.e1[2]
        shearStrain = s.e1[3]
        shearStress = s.s1[3]
        eplEquiv = s.getInternalStateVariableValue("EquivalentPlasticStrain")
        porosity = s.getInternalStateVariableValue("Porosity")
        volRatio = s.getInternalStateVariableValue("VolumeRatio")
        pc = s.getInternalStateVariableValue("PreConsolidationPressure")
        eplV = s.getInternalStateVariableValue("PlasticVolumetricStrain")

        results[0][k][i + 1] = vMstress  # / pc
        results[1][k][i + 1] = pressure  # / pc,  + pamb
        results[2][k][i + 1] = volRatio
        results[3][k][i + 1] = porosity
        results[4][k][i + 1] = shearStrain
        results[5][k][i + 1] = shearStress
        results[6][k][i + 1] = pc
        results[7][k][i + 1] = eplV

    # save the (asymptotic) values of q,p -> should be on the CSL
    qEndvalues[k] = results[0][k][nTime - 1]
    pEndvalues[k] = results[1][k][nTime - 1]
    # calculate the final yield surfaces in the p-q-space
    pRangeFinal[k] = np.linspace(-0, pc, nPoints)
    qFunctFinal[k] = M * np.sqrt((pc - pRangeFinal[k]) * pRangeFinal[k])
    # calculate the theoretical CSL
    qFunctCSL = M * pRangeFinal[k]

# plots
fig, ax = plt.subplots()
ax.set_title("Total volume/solid volume over time")
for k in range(runs):
    ax.plot(
        ltime, results[2][k], label=prelabel + "%.2f" % (valueList[k])
    )  # OCR: pc0/valueList[k]
ax.set_xlabel("$t$ / s")
ax.set_ylabel("volume ratio")
# ax.set_ylim(1.7875,1.788)
ax.grid()
ax.legend()
fig.savefig("out/ModCamClay_ParamStudy_VolumeRatio.pdf")

fig, ax = plt.subplots()
ax.set_title("Porosity over time")
for k in range(runs):
    ax.plot(ltime, results[3][k], label=prelabel + "%.2f" % (valueList[k]))
ax.set_xlabel("$t$ / s")
ax.set_ylabel("$\phi$")
ax.grid()
ax.legend()
fig.savefig("out/ModCamClay_ParamStudy_Porosity.pdf")

fig, ax = plt.subplots()
# ax.set_title('Shear stress over shear strain')
for k in range(runs):
    ax.plot(results[4][k], results[5][k], label=prelabel + "%.2f" % (valueList[k]))
ax.set_xlabel("$\epsilon_{xy}$")
ax.set_ylabel("$\sigma_{xy}$ / MPa")
ax.grid()
ax.legend()
fig.savefig("out/ModCamClay_ParamStudy_ShearCurves.pdf")

fig, ax = plt.subplots()
# ax.set_title('Plastic volumetric strain over shear strain')
for k in range(runs):
    ax.plot(
        results[4][k], results[7][k] * 100, label=prelabel + "%.2f" % (valueList[k])
    )
ax.set_xlabel("$\epsilon_{xy}$")
ax.set_ylabel("${\epsilon}_p^V$ / %")
ax.grid()
ax.legend()
fig.savefig("out/ModCamClay_ParamStudy_eplVCurves.pdf")

fig, ax = plt.subplots()
ax.set_title("Hydrostatic pressure and von-Mises stress over time")
for k in range(runs):
    ax.plot(ltime, results[0][k], label=prelabel + "%.2f" % (valueList[k]))
    ax.plot(
        ltime,
        results[1][k],
        label=prelabel + "%.2f" % (valueList[k]),
        linestyle="dashed",
    )
ax.set_xlabel("$t$ / s")
ax.set_ylabel("$q, p$ / MPa")
ax.grid()
fig.savefig("out/ModCamClay_ParamStudy_qpCurves.pdf")

fig, ax = plt.subplots()
ax.set_title("Hydrostatic pressure and pre-consolidation pressure over time")
for k in range(runs):
    ax.plot(
        ltime,
        results[1][k],
        label=prelabel + "%.2f" % (valueList[k]),
        linestyle="dashed",
    )
    ax.plot(ltime, results[6][k], label=prelabel + "%.2f" % (valueList[k]))
ax.set_xlabel("$t$ / s")
ax.set_ylabel("$p, p_c$ / MPa")
ax.grid()
ax.legend(loc="upper left")
fig.savefig("out/ModCamClay_ParamStudy_PressCurves.pdf")

fig, ax = plt.subplots()
# ax.set_title('Values of p and q at the initial and final yield surfaces')
ax.plot(pRange, qFunct, color="grey", label="Initial")
ax.plot(pRange, M * pRange, color="black", label="CSL")
for k in range(runs):
    # final yield surface
    ax.scatter(pRangeFinal[k], qFunctFinal[k], label=prelabel + "%.2f" % (valueList[k]))
    # stress trajectory
    ax.plot(results[1][k], results[0][k], label=prelabel + "%.2f" % (valueList[k]))
ax.set_xlabel("$p$ / MPa")
ax.set_ylabel("$q$ / MPa")
ax.grid()
# legend outside by resizing the box and putting the legend relative to that
chartBox = ax.get_position()
ax.set_position([chartBox.x0, chartBox.y0, chartBox.width * 0.85, chartBox.height])
ax.legend(loc="upper center", bbox_to_anchor=(1.15, 0.80))
# ax.legend()
# ax.legend(loc='upper left')
fig.savefig("out/ModCamClay_ParamStudy_YieldSurface.pdf")

if variedParameter == "p":
    fig, ax = plt.subplots()
    ax.set_title("Values of p and q at simulation end - Critical State Line")
    ax.plot(pEndvalues, M * pEndvalues, label="analytical CSL")
    ax.plot(pEndvalues, qEndvalues, label="calculated CSL", zorder=1)
    ax.scatter(pEndvalues, qEndvalues, color="r", zorder=2)
    ax.set_xlabel("$p$ / MPa")
    ax.set_ylabel("$q$ / MPa")
    ax.grid()
    ax.legend()
    fig.savefig("out/ModCamClay_ParamStudy_CSL.pdf")


plt.show()
