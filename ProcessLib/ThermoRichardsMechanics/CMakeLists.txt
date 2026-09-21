set(_directories . ConstitutiveCommon ConstitutiveStress_StrainTemperature)
if(OGS_USE_MFRONT)
    list(APPEND _directories
         ConstitutiveStressSaturation_StrainPressureTemperature
    )
endif()

if(WIN32)
    # Exceeds library limit of 65535 objects
    set(_static STATIC)
endif()

ogs_add_process(ThermoRichardsMechanics ${_static} DIRECTORIES ${_directories})
