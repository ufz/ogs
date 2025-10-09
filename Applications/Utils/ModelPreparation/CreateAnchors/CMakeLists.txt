set(EXE_SOURCES CreateAnchors.cpp ComputeIntersections.cpp)

ogs_add_executable(CreateAnchors ${EXE_SOURCES})

target_link_libraries(CreateAnchors tclap nlohmann_json::nlohmann_json VTK::FiltersGeneral VTK::FiltersSources NumLib)
string(REPLACE "." "" _vtk_version_stripped "${VTK_VERSION}")
target_compile_definitions(CreateAnchors PRIVATE VTK_VERSION_STRIPPED=${_vtk_version_stripped})

install(TARGETS CreateAnchors RUNTIME DESTINATION bin)
