/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef ASSEMBLERLIB_GLOBALSETUP_H_
#define ASSEMBLERLIB_GLOBALSETUP_H_

#include <functional>

namespace AssemblerLib
{

/// The GlobalSetup collects vector and matrix builder and corresponding global
/// loop executor.
template <typename VectorMatrixBuilder, typename Executor, typename LinearSolver_>
struct GlobalSetup
{
    typedef typename VectorMatrixBuilder::VectorType VectorType;
    typedef typename VectorMatrixBuilder::MatrixType MatrixType;

    using LinearSolver = LinearSolver_;

    template <typename... Args>
    static
    VectorType* createVector(Args&& ... args)
    {
        return VectorMatrixBuilder::createVector(std::forward<Args>(args)...);
    }

    template <typename... Args>
    static
    MatrixType* createMatrix(Args&& ... args)
    {
        return VectorMatrixBuilder::createMatrix(std::forward<Args>(args)...);
    }

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
    void transform(Args&& ... args)
    {
        return Executor::transform(std::forward<Args>(args)...);
    }

    //! Do not create any instances; this struct only has static members.
    GlobalSetup() = delete;
};

}   // namespace AssemblerLib

#endif  // ASSEMBLERLIB_GLOBALSETUP_H_
