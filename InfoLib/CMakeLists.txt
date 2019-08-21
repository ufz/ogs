foreach(lib Git CMake Test)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${lib}Info.cpp.in
                   ${CMAKE_CURRENT_BINARY_DIR}/${lib}Info.cpp @ONLY)

    add_library(${lib}InfoLib ${CMAKE_CURRENT_BINARY_DIR}/${lib}Info.cpp
                              ${lib}Info.h)

    include(GenerateExportHeader)
    generate_export_header(${lib}InfoLib)
    target_include_directories(${lib}InfoLib PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
endforeach(lib)
