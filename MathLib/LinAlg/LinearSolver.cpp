/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "LinearSolver.h"

namespace MathLib
{

template <typename Solver>
std::unique_ptr<Solver> createLinearSolver(
    BaseLib::ConfigTree const* const config)
{
    return std::unique_ptr<Solver>{new Solver{"", config}};
}

} // MathLib


// does explicit template instantiation of createLinearSolver<>()
#define SPECIALIZE_CREATE_LINEAR_SOLVER(SLV) \
    namespace MathLib { \
    template std::unique_ptr<SLV> \
    createLinearSolver<SLV>(BaseLib::ConfigTree const*const); \
    }


#ifdef OGS_USE_LIS

#include "MathLib/LinAlg/EigenLis/EigenLisLinearSolver.h"

SPECIALIZE_CREATE_LINEAR_SOLVER(MathLib::EigenLisLinearSolver)

#elif defined(USE_PETSC)

#include "MathLib/LinAlg/PETSc/PETScLinearSolver.h"

SPECIALIZE_CREATE_LINEAR_SOLVER(MathLib::PETScLinearSolver)

#elif defined(OGS_USE_EIGEN)

#include "MathLib/LinAlg/Eigen/EigenLinearSolver.h"

SPECIALIZE_CREATE_LINEAR_SOLVER(MathLib::EigenLinearSolver)

#endif
