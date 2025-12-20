// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

// TODO merge that with MatrixVectorTraits?

#pragma once

#include <Eigen/Core>
#include <initializer_list>

#include "MatrixVectorTraits.h"

#ifdef USE_PETSC

// Global PETScMatrix/PETScVector //////////////////////////////////////////
namespace MathLib
{

class PETScVector;
class PETScMatrix;

void setVector(PETScVector& v, std::initializer_list<double> values);

void setVector(PETScVector& v,
               MatrixVectorTraits<PETScVector>::Index const index,
               double const value);

void setMatrix(PETScMatrix& m, Eigen::MatrixXd const& tmp);

void addToMatrix(PETScMatrix& m, std::initializer_list<double> values);

void setMatrix(PETScMatrix& m, std::initializer_list<double> values);

}  // namespace MathLib

#else

// Sparse global EigenMatrix/EigenVector
// //////////////////////////////////////////

namespace MathLib
{

class EigenVector;
class EigenMatrix;

void setVector(EigenVector& v, std::initializer_list<double> values);

void setVector(EigenVector& v,
               MatrixVectorTraits<EigenVector>::Index const index,
               double const value);

void setMatrix(EigenMatrix& m, std::initializer_list<double> values);

void setMatrix(EigenMatrix& m, Eigen::MatrixXd const& tmp);

void addToMatrix(EigenMatrix& m, std::initializer_list<double> values);

}  // namespace MathLib

#endif
