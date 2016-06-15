/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef NUMLIB_GLOBAL_SETUP_H_
#define NUMLIB_GLOBAL_SETUP_H_

#include <functional>

#include "MathLib/LinAlg/GlobalMatrixVectorTypes.h"

namespace NumLib
{

/// The GlobalSetup collects vector and matrix builder and corresponding global
/// loop executor.
//template <typename VectorMatrixBuilder, typename Executor, typename LinearSolver_>
template <typename Executor>
struct GlobalSetup
{
    using VectorType = detail::GlobalVectorType;
    using MatrixType = detail::GlobalMatrixType;

    using LinearSolver = detail::LinearSolverType;

//    template <typename... Args>
//    static
//    VectorType* createVector(Args&& ... args)
//    {
//        return VectorMatrixBuilder::createVector(std::forward<Args>(args)...);
//    }

//    template <typename... Args>
//    static
//    MatrixType* createMatrix(Args&& ... args)
//    {
//        return VectorMatrixBuilder::createMatrix(std::forward<Args>(args)...);
//    }

    template <typename... Args>
    static
    void executeDereferenced(Args&& ... args)
    {
        return Executor::executeDereferenced(std::forward<Args>(args)...);
    }

    template <typename... Args>
    static
    void executeMemberDereferenced(Args&& ... args)
    {
        return Executor::executeMemberDereferenced(std::forward<Args>(args)...);
    }

    template <typename... Args>
    static
    void transformDereferenced(Args&& ... args)
    {
        return Executor::transformDereferenced(std::forward<Args>(args)...);
    }

    //! Do not create any instances; this struct only has static members.
    GlobalSetup() = delete;
};

}   // namespace NumLib

#endif  // NUMLIB_GLOBAL_SETUP_H_
