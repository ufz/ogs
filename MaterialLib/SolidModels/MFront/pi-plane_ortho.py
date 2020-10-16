from math import pi, cos, sin, sqrt
from tfel.material import projectOnPiPlane
import mtest

divisions = 1000
for theta in [
    pi * (-1.0 + 2.0 * float(i) / (float(divisions) - 1.0)) for i in range(0, divisions)
]:
    # for theta in [-1.3010636242139548]:
    em = 5.0e-3
    npas = 100
    tmax = 1
    c = cos(theta)
    s = sin(theta)
    m = mtest.MTest()
    mtest.setVerboseMode(mtest.VerboseLevel.VERBOSE_QUIET)
    m.setMaximumNumberOfSubSteps(10)
    m.setBehaviour("generic", "src/libBehaviour.so", "MohrCoulombAbboSloanOrtho")
    m.setExternalStateVariable("Temperature", 293.15)
    m.setImposedStrain("EXX", {0: 0, tmax: em * c})
    m.setImposedStrain("EYY", {0: 0, tmax: em * s})
    m.setNonLinearConstraint("SXX+SYY+SZZ", "Stress")
    m.setNonLinearConstraint("SXY", "Stress")
    m.setNonLinearConstraint("SXZ", "Stress")
    m.setNonLinearConstraint("SYZ", "Stress")
    m.setMaterialProperty("YoungModulus1", 150.0e3)
    m.setMaterialProperty("YoungModulus2", 50.0e3)
    m.setMaterialProperty("YoungModulus3", 250.0e3)
    m.setMaterialProperty("PoissonRatio12", 0.13)
    m.setMaterialProperty("PoissonRatio23", 0.24)
    m.setMaterialProperty("PoissonRatio13", 0.18)
    m.setMaterialProperty("ShearModulus12", 60e3)
    m.setMaterialProperty("ShearModulus23", 30e3)
    m.setMaterialProperty("ShearModulus13", 180e3)
    m.setMaterialProperty("Cohesion", 3.0e1)
    m.setMaterialProperty("FrictionAngle", 30.0)
    m.setMaterialProperty("DilatancyAngle", 10.0)
    m.setMaterialProperty("TransitionAngle", 29.0)
    m.setMaterialProperty("TensionCutOffParameter", 1.0e1)
    s = mtest.MTestCurrentState()
    wk = mtest.MTestWorkSpace()
    m.completeInitialisation()
    m.initializeCurrentState(s)
    m.initializeWorkSpace(wk)
    ltime = [float((tmax / (npas - 1))) * i for i in range(npas)]
    plas = 0
    plas_tol = 1e-10
    for i in range(npas - 1):
        m.execute(s, wk, ltime[i], ltime[i + 1])
        p = s.getInternalStateVariableValue("EquivalentPlasticStrain")
        s0, s1 = projectOnPiPlane(s.s1[0], s.s1[1], s.s1[2])
        if p > plas_tol:
            print("0. 0. " + str(s0) + " " + str(s1))
            plas += 1
            if plas > 1:
                break
        else:
            print(str(s0) + " " + str(s1) + " 0. 0.")
