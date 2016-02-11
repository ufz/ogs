#pragma once

#include <Eigen/SparseCore>

using IndexType = int;
using Matrix = Eigen::SparseMatrix<double, 0, IndexType>;
using Vector = Eigen::VectorXd;


enum class NonlinearSolverTag : bool { Picard, Newton };

enum class ODESystemTag : char
{
    FirstOrderImplicitQuasilinear
};
