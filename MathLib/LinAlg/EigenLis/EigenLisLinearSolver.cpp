/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "EigenLisLinearSolver.h"

#ifdef _OPENMP
#include <omp.h>
#endif
#include <logog/include/logog.hpp>

#include "BaseLib/ConfigTreeNew.h"
#include "MathLib/LinAlg/Eigen/EigenMatrix.h"
#include "MathLib/LinAlg/Eigen/EigenVector.h"
#include "MathLib/LinAlg/Lis/LisMatrix.h"
#include "MathLib/LinAlg/Lis/LisVector.h"
#include "MathLib/LinAlg/Lis/LisLinearSolver.h"

namespace MathLib
{
EigenLisLinearSolver::EigenLisLinearSolver(
    EigenMatrix& A,
    const std::string /*solver_name*/,
    BaseLib::ConfigTree const* const option)
    : _A(A), _lis_option(option)
{
}

void EigenLisLinearSolver::solve(EigenVector &b_, EigenVector &x_)
{
    static_assert(EigenMatrix::RawMatrixType::IsRowMajor,
                  "Sparse matrix is required to be in row major storage.");
    auto &A = _A.getRawMatrix();
    auto &b = b_.getRawVector();
    auto &x = x_.getRawVector();

    if (!A.isCompressed())
        A.makeCompressed();
    int nnz = A.nonZeros();
    int* ptr = A.outerIndexPtr();
    int* col = A.innerIndexPtr();
    double* data = A.valuePtr();
    LisMatrix lisA(_A.getNRows(), nnz, ptr, col, data);
    LisVector lisb(b.rows(), b.data());
    LisVector lisx(x.rows(), x.data());

    LisLinearSolver lissol(lisA);
    lissol.setOption(_lis_option);
    lissol.solve(lisb, lisx);

    for (std::size_t i=0; i<lisx.size(); i++)
        x[i] = lisx[i];
}

} //MathLib
