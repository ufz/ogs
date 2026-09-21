ogs_add_process(ThermoHydroMechanics GENERATE_EXPORT_HEADER)

# Disable Windows DLL export all symbols for this target to reduce symbol count
# and avoid LNK1189 error (65535 object limit exceeded)
set_target_properties(
    ThermoHydroMechanics PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS OFF
)

# Set all symbols to hidden visibility by default
set_target_properties(
    ThermoHydroMechanics
    PROPERTIES CXX_VISIBILITY_PRESET hidden C_VISIBILITY_PRESET hidden
               VISIBILITY_INLINES_HIDDEN ON
)
