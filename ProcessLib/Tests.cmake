# This contains multi-process notebooks
if(TARGET ThermoHydroMechanics
   AND TARGET ThermoRichardsMechanics
   AND TARGET TH2M
   AND NOT OGS_USE_PETSC
)
    NotebookTest(
        NOTEBOOKFILE
        ThermoHydroMechanics/Linear/Point_injection/SaturatedPointheatsource.ipynb
        RUNTIME 1800
        PROPERTIES PROCESSORS 4
    )
endif()
