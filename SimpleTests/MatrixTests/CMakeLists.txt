# Create the executable
add_executable(DenseGaussEliminationChecker
    DenseGaussEliminationChecker.cpp
    ${SOURCES}
    ${HEADERS}
)
set_target_properties(DenseGaussEliminationChecker PROPERTIES FOLDER SimpleTests)
target_link_libraries(DenseGaussEliminationChecker
    logog
    BaseLib
    MathLib
)
