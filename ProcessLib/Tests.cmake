# This contains multi-process notebooks
if(TARGET ThermoHydroMechanics
   AND TARGET ThermoRichardsMechanics
   AND TARGET TH2M
   AND NOT OGS_USE_PETSC
)
    NotebookTest(
        NOTEBOOKFILE
        ThermoHydroMechanics/Linear/Point_injection/SaturatedPointheatsource.ipynb
        RUNTIME
        1800
    )
    if(TEST
       nb-ThermoHydroMechanics/Linear/Point_injection/SaturatedPointheatsource-LARGE
    )
        set_tests_properties(
            nb-ThermoHydroMechanics/Linear/Point_injection/SaturatedPointheatsource-LARGE
            PROPERTIES PROCESSORS 4
        )
    endif()
endif()
