add_executable(xdmfdiff xdmfdiff.cpp)
target_link_libraries(xdmfdiff OgsXdmf std::filesystem tclap)
install(PROGRAMS $<TARGET_FILE:xdmfdiff> DESTINATION bin)
