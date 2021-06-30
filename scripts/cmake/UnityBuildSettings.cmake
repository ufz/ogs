if(NOT OGS_USE_UNITY_BUILDS OR ${CMAKE_VERSION} VERSION_LESS 3.16)
    return()
endif()

set_target_properties(BaseLib PROPERTIES UNITY_BUILD_BATCH_SIZE 8)
set_target_properties(GeoLib PROPERTIES UNITY_BUILD_BATCH_SIZE 40)
if(TARGET MaterialLib)
    set_target_properties(MaterialLib PROPERTIES UNITY_BUILD_BATCH_SIZE 20)
endif()
set_target_properties(MathLib PROPERTIES UNITY_BUILD_BATCH_SIZE 10)
set_target_properties(MeshLib PROPERTIES UNITY_BUILD_BATCH_SIZE 20)
# set_target_properties(ProcessLib PROPERTIES UNITY_BUILD_BATCH_SIZE 80) #
# breaks!

if(TARGET testrunner)
    set_target_properties(testrunner PROPERTIES UNITY_BUILD ON)
endif()
