import mtest

em = 1.0e-2
npas = 100
tmax = 1 * 6
m = mtest.MTest()
mtest.setVerboseMode(mtest.VerboseLevel.VERBOSE_QUIET)
m.setMaximumNumberOfSubSteps(10)
m.setBehaviour("generic", "src/libBehaviour.so", "MohrCoulombAbboSloanOrtho")
m.setExternalStateVariable("Temperature", 293.15)
m.setImposedStrain("EXX", {0: 0, tmax / 6: em, tmax / 3: 0})
m.setImposedStrain("EYY", {tmax / 3: 0, tmax / 2: em, 2 * tmax / 3: 0})
m.setImposedStrain("EZZ", {2 * tmax / 3: 0, 5 * tmax / 6: em, tmax: 0})
m.setMaterialProperty("YoungModulus1", 150.0e3)
m.setMaterialProperty("YoungModulus2", 50.0e3)
m.setMaterialProperty("YoungModulus3", 250.0e3)
m.setMaterialProperty("PoissonRatio12", 0.13)
m.setMaterialProperty("PoissonRatio23", 0.24)
m.setMaterialProperty("PoissonRatio13", 0.18)
m.setMaterialProperty("ShearModulus12", 60e3)
m.setMaterialProperty("ShearModulus23", 30e3)
m.setMaterialProperty("ShearModulus13", 180e3)
m.setMaterialProperty("Cohesion", 3.0e19)
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
for i in range(npas - 1):
    m.execute(s, wk, ltime[i], ltime[i + 1])
    print(
        str(s.e1[0])
        + " "
        + str(s.e1[1])
        + " "
        + str(s.e1[2])
        + " "
        + str(s.e1[3])
        + " "
        + str(s.e1[4])
        + " "
        + str(s.e1[5])
        + " "
        + str(s.s1[0])
        + " "
        + str(s.s1[1])
        + " "
        + str(s.s1[2])
        + " "
        + str(s.s1[3])
        + " "
        + str(s.s1[4])
        + " "
        + str(s.s1[5])
    )
