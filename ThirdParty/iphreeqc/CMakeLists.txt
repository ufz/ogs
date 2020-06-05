project(IPhreeqc CXX C)
set(IPHREEQC_VERSION 3.5.0)

set(SOURCES
    src/src/CSelectedOutput.cpp
    src/src/fwrap.cpp
    src/src/IPhreeqc.cpp
    src/src/IPhreeqc_interface_F.cpp
    src/src/IPhreeqcLib.cpp
    src/src/Var.c
    src/src/phreeqcpp/advection.cpp
    src/src/phreeqcpp/basicsubs.cpp
    src/src/phreeqcpp/cl1.cpp
    src/src/phreeqcpp/cvdense.cpp
    src/src/phreeqcpp/cvode.cpp
    src/src/phreeqcpp/cxxKinetics.cxx
    src/src/phreeqcpp/cxxMix.cxx
    src/src/phreeqcpp/dense.cpp
    src/src/phreeqcpp/Dictionary.cpp
    src/src/phreeqcpp/dumper.cpp
    src/src/phreeqcpp/Exchange.cxx
    src/src/phreeqcpp/ExchComp.cxx
    src/src/phreeqcpp/GasComp.cxx
    src/src/phreeqcpp/gases.cpp
    src/src/phreeqcpp/GasPhase.cxx
    src/src/phreeqcpp/input.cpp
    src/src/phreeqcpp/integrate.cpp
    src/src/phreeqcpp/inverse.cpp
    src/src/phreeqcpp/ISolution.cxx
    src/src/phreeqcpp/ISolutionComp.cxx
    src/src/phreeqcpp/isotopes.cpp
    src/src/phreeqcpp/kinetics.cpp
    src/src/phreeqcpp/KineticsComp.cxx
    src/src/phreeqcpp/mainsubs.cpp
    src/src/phreeqcpp/model.cpp
    src/src/phreeqcpp/NameDouble.cxx
    src/src/phreeqcpp/NumKeyword.cxx
    src/src/phreeqcpp/nvector.cpp
    src/src/phreeqcpp/nvector_serial.cpp
    src/src/phreeqcpp/parse.cpp
    src/src/phreeqcpp/PBasic.cpp
    src/src/phreeqcpp/phqalloc.cpp
    src/src/phreeqcpp/Phreeqc.cpp
    src/src/phreeqcpp/PHRQ_io_output.cpp
    src/src/phreeqcpp/pitzer.cpp
    src/src/phreeqcpp/pitzer_structures.cpp
    src/src/phreeqcpp/PPassemblage.cxx
    src/src/phreeqcpp/PPassemblageComp.cxx
    src/src/phreeqcpp/prep.cpp
    src/src/phreeqcpp/Pressure.cxx
    src/src/phreeqcpp/print.cpp
    src/src/phreeqcpp/Reaction.cxx
    src/src/phreeqcpp/read.cpp
    src/src/phreeqcpp/ReadClass.cxx
    src/src/phreeqcpp/readtr.cpp
    src/src/phreeqcpp/runner.cpp
    src/src/phreeqcpp/SelectedOutput.cpp
    src/src/phreeqcpp/Serializer.cxx
    src/src/phreeqcpp/sit.cpp
    src/src/phreeqcpp/smalldense.cpp
    src/src/phreeqcpp/Solution.cxx
    src/src/phreeqcpp/SolutionIsotope.cxx
    src/src/phreeqcpp/spread.cpp
    src/src/phreeqcpp/SS.cxx
    src/src/phreeqcpp/SSassemblage.cxx
    src/src/phreeqcpp/SScomp.cxx
    src/src/phreeqcpp/step.cpp
    src/src/phreeqcpp/StorageBin.cxx
    src/src/phreeqcpp/StorageBinList.cpp
    src/src/phreeqcpp/structures.cpp
    src/src/phreeqcpp/sundialsmath.cpp
    src/src/phreeqcpp/Surface.cxx
    src/src/phreeqcpp/SurfaceCharge.cxx
    src/src/phreeqcpp/SurfaceComp.cxx
    src/src/phreeqcpp/System.cxx
    src/src/phreeqcpp/tally.cpp
    src/src/phreeqcpp/Temperature.cxx
    src/src/phreeqcpp/tidy.cpp
    src/src/phreeqcpp/transport.cpp
    src/src/phreeqcpp/Use.cpp
    src/src/phreeqcpp/UserPunch.cpp
    src/src/phreeqcpp/utilities.cpp
    src/src/phreeqcpp/common/Parser.cxx
    src/src/phreeqcpp/common/PHRQ_base.cxx
    src/src/phreeqcpp/common/PHRQ_io.cpp
    src/src/phreeqcpp/common/Utils.cxx
    src/src/phreeqcpp/PhreeqcKeywords/Keywords.cpp
)

# compile Var.c as c++
SET_SOURCE_FILES_PROPERTIES(src/src/Var.c PROPERTIES LANGUAGE CXX)

add_library(iphreeqc STATIC ${SOURCES})
target_include_directories(iphreeqc BEFORE
    PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/src/src/phreeqcpp/common
    ${CMAKE_CURRENT_SOURCE_DIR}/src/src/phreeqcpp/PhreeqcKeywords
    PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src/src
    ${CMAKE_CURRENT_SOURCE_DIR}/src/src/phreeqcpp
)
# Exclude iphreeqc target from clang-tidy tests because it handles the above
# mentioned 'src/src/Var.c' file as c, not c++.
set_target_properties(iphreeqc PROPERTIES CXX_CLANG_TIDY "")
