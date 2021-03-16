set(SOURCES
    ${iphreeqc_SOURCE_DIR}/src/CSelectedOutput.cpp
    ${iphreeqc_SOURCE_DIR}/src/fwrap.cpp
    ${iphreeqc_SOURCE_DIR}/src/IPhreeqc.cpp
    ${iphreeqc_SOURCE_DIR}/src/IPhreeqc_interface_F.cpp
    ${iphreeqc_SOURCE_DIR}/src/IPhreeqcLib.cpp
    ${iphreeqc_SOURCE_DIR}/src/Var.c
    ${iphreeqc_SOURCE_DIR}/src/phreeqcpp/advection.cpp
    ${iphreeqc_SOURCE_DIR}/src/phreeqcpp/basicsubs.cpp
    ${iphreeqc_SOURCE_DIR}/src/phreeqcpp/cl1.cpp
    ${iphreeqc_SOURCE_DIR}/src/phreeqcpp/cvdense.cpp
    ${iphreeqc_SOURCE_DIR}/src/phreeqcpp/cvode.cpp
    ${iphreeqc_SOURCE_DIR}/src/phreeqcpp/cxxKinetics.cxx
    ${iphreeqc_SOURCE_DIR}/src/phreeqcpp/cxxMix.cxx
    ${iphreeqc_SOURCE_DIR}/src/phreeqcpp/dense.cpp
    ${iphreeqc_SOURCE_DIR}/src/phreeqcpp/Dictionary.cpp
    ${iphreeqc_SOURCE_DIR}/src/phreeqcpp/dumper.cpp
    ${iphreeqc_SOURCE_DIR}/src/phreeqcpp/Exchange.cxx
    ${iphreeqc_SOURCE_DIR}/src/phreeqcpp/ExchComp.cxx
    ${iphreeqc_SOURCE_DIR}/src/phreeqcpp/GasComp.cxx
    ${iphreeqc_SOURCE_DIR}/src/phreeqcpp/gases.cpp
    ${iphreeqc_SOURCE_DIR}/src/phreeqcpp/GasPhase.cxx
    ${iphreeqc_SOURCE_DIR}/src/phreeqcpp/input.cpp
    ${iphreeqc_SOURCE_DIR}/src/phreeqcpp/integrate.cpp
    ${iphreeqc_SOURCE_DIR}/src/phreeqcpp/inverse.cpp
    ${iphreeqc_SOURCE_DIR}/src/phreeqcpp/ISolution.cxx
    ${iphreeqc_SOURCE_DIR}/src/phreeqcpp/ISolutionComp.cxx
    ${iphreeqc_SOURCE_DIR}/src/phreeqcpp/isotopes.cpp
    ${iphreeqc_SOURCE_DIR}/src/phreeqcpp/kinetics.cpp
    ${iphreeqc_SOURCE_DIR}/src/phreeqcpp/KineticsComp.cxx
    ${iphreeqc_SOURCE_DIR}/src/phreeqcpp/mainsubs.cpp
    ${iphreeqc_SOURCE_DIR}/src/phreeqcpp/model.cpp
    ${iphreeqc_SOURCE_DIR}/src/phreeqcpp/NameDouble.cxx
    ${iphreeqc_SOURCE_DIR}/src/phreeqcpp/NumKeyword.cxx
    ${iphreeqc_SOURCE_DIR}/src/phreeqcpp/nvector.cpp
    ${iphreeqc_SOURCE_DIR}/src/phreeqcpp/nvector_serial.cpp
    ${iphreeqc_SOURCE_DIR}/src/phreeqcpp/parse.cpp
    ${iphreeqc_SOURCE_DIR}/src/phreeqcpp/PBasic.cpp
    ${iphreeqc_SOURCE_DIR}/src/phreeqcpp/phqalloc.cpp
    ${iphreeqc_SOURCE_DIR}/src/phreeqcpp/Phreeqc.cpp
    ${iphreeqc_SOURCE_DIR}/src/phreeqcpp/PHRQ_io_output.cpp
    ${iphreeqc_SOURCE_DIR}/src/phreeqcpp/pitzer.cpp
    ${iphreeqc_SOURCE_DIR}/src/phreeqcpp/pitzer_structures.cpp
    ${iphreeqc_SOURCE_DIR}/src/phreeqcpp/PPassemblage.cxx
    ${iphreeqc_SOURCE_DIR}/src/phreeqcpp/PPassemblageComp.cxx
    ${iphreeqc_SOURCE_DIR}/src/phreeqcpp/prep.cpp
    ${iphreeqc_SOURCE_DIR}/src/phreeqcpp/Pressure.cxx
    ${iphreeqc_SOURCE_DIR}/src/phreeqcpp/print.cpp
    ${iphreeqc_SOURCE_DIR}/src/phreeqcpp/Reaction.cxx
    ${iphreeqc_SOURCE_DIR}/src/phreeqcpp/read.cpp
    ${iphreeqc_SOURCE_DIR}/src/phreeqcpp/ReadClass.cxx
    ${iphreeqc_SOURCE_DIR}/src/phreeqcpp/readtr.cpp
    ${iphreeqc_SOURCE_DIR}/src/phreeqcpp/runner.cpp
    ${iphreeqc_SOURCE_DIR}/src/phreeqcpp/SelectedOutput.cpp
    ${iphreeqc_SOURCE_DIR}/src/phreeqcpp/Serializer.cxx
    ${iphreeqc_SOURCE_DIR}/src/phreeqcpp/sit.cpp
    ${iphreeqc_SOURCE_DIR}/src/phreeqcpp/smalldense.cpp
    ${iphreeqc_SOURCE_DIR}/src/phreeqcpp/Solution.cxx
    ${iphreeqc_SOURCE_DIR}/src/phreeqcpp/SolutionIsotope.cxx
    ${iphreeqc_SOURCE_DIR}/src/phreeqcpp/spread.cpp
    ${iphreeqc_SOURCE_DIR}/src/phreeqcpp/SS.cxx
    ${iphreeqc_SOURCE_DIR}/src/phreeqcpp/SSassemblage.cxx
    ${iphreeqc_SOURCE_DIR}/src/phreeqcpp/SScomp.cxx
    ${iphreeqc_SOURCE_DIR}/src/phreeqcpp/step.cpp
    ${iphreeqc_SOURCE_DIR}/src/phreeqcpp/StorageBin.cxx
    ${iphreeqc_SOURCE_DIR}/src/phreeqcpp/StorageBinList.cpp
    ${iphreeqc_SOURCE_DIR}/src/phreeqcpp/structures.cpp
    ${iphreeqc_SOURCE_DIR}/src/phreeqcpp/sundialsmath.cpp
    ${iphreeqc_SOURCE_DIR}/src/phreeqcpp/Surface.cxx
    ${iphreeqc_SOURCE_DIR}/src/phreeqcpp/SurfaceCharge.cxx
    ${iphreeqc_SOURCE_DIR}/src/phreeqcpp/SurfaceComp.cxx
    ${iphreeqc_SOURCE_DIR}/src/phreeqcpp/System.cxx
    ${iphreeqc_SOURCE_DIR}/src/phreeqcpp/tally.cpp
    ${iphreeqc_SOURCE_DIR}/src/phreeqcpp/Temperature.cxx
    ${iphreeqc_SOURCE_DIR}/src/phreeqcpp/tidy.cpp
    ${iphreeqc_SOURCE_DIR}/src/phreeqcpp/transport.cpp
    ${iphreeqc_SOURCE_DIR}/src/phreeqcpp/Use.cpp
    ${iphreeqc_SOURCE_DIR}/src/phreeqcpp/UserPunch.cpp
    ${iphreeqc_SOURCE_DIR}/src/phreeqcpp/utilities.cpp
    ${iphreeqc_SOURCE_DIR}/src/phreeqcpp/common/Parser.cxx
    ${iphreeqc_SOURCE_DIR}/src/phreeqcpp/common/PHRQ_base.cxx
    ${iphreeqc_SOURCE_DIR}/src/phreeqcpp/common/PHRQ_io.cpp
    ${iphreeqc_SOURCE_DIR}/src/phreeqcpp/common/Utils.cxx
    ${iphreeqc_SOURCE_DIR}/src/phreeqcpp/PhreeqcKeywords/Keywords.cpp
)

# compile Var.c as c++
SET_SOURCE_FILES_PROPERTIES(
    ${iphreeqc_SOURCE_DIR}/src/Var.c PROPERTIES LANGUAGE CXX
)

add_library(iphreeqc STATIC ${SOURCES})
target_include_directories(
    iphreeqc
    PUBLIC ${iphreeqc_SOURCE_DIR}/src/phreeqcpp/common
           ${iphreeqc_SOURCE_DIR}/src/phreeqcpp/PhreeqcKeywords
           ${iphreeqc_SOURCE_DIR}/src/phreeqcpp ${iphreeqc_SOURCE_DIR}/src
)
target_compile_definitions(iphreeqc PUBLIC LDBLE=double)
# Exclude iphreeqc target from clang-tidy tests because it handles the above
# mentioned 'src/src/Var.c' file as c, not c++.
set_target_properties(iphreeqc PROPERTIES CXX_CLANG_TIDY "")
