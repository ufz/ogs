/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

// TODO merge that with MatrixVectorTraits?

#pragma once

#include <initializer_list>
#include <Eigen/Core>
#include "MatrixVectorTraits.h"

#ifdef USE_PETSC

// Global PETScMatrix/PETScVector //////////////////////////////////////////
namespace MathLib
{

class PETScVector;
class PETScMatrix;

void setVector(PETScVector& v,
               std::initializer_list<double> values);

void setVector(PETScVector& v, MatrixVectorTraits<PETScVector>::Index const index,
               double const value);

void setMatrix(PETScMatrix& m, Eigen::MatrixXd const& tmp);

void addToMatrix(PETScMatrix& m,
                 std::initializer_list<double> values);

void setMatrix(PETScMatrix& m,
               std::initializer_list<double> values);

} // namespace MathLib


#elif defined(OGS_USE_EIGEN)

// Sparse global EigenMatrix/EigenVector //////////////////////////////////////////

namespace MathLib
{

class EigenVector;
class EigenMatrix;

void setVector(EigenVector& v,
               std::initializer_list<double> values);

void setVector(EigenVector& v, MatrixVectorTraits<EigenVector>::Index const index,
               double const value);

void setMatrix(EigenMatrix& m,
               std::initializer_list<double> values);

void setMatrix(EigenMatrix& m, Eigen::MatrixXd const& tmp);

void addToMatrix(EigenMatrix& m,
                 std::initializer_list<double> values);

} // namespace MathLib

#endif // OGS_USE_EIGEN
