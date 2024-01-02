/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <Eigen/Core>

namespace MathLib
{
template <int D, typename M>
class EigenBlockMatrixViewFunctor
{
public:
    static constexpr int rows = M::RowsAtCompileTime == Eigen::Dynamic
                                    ? Eigen::Dynamic
                                    : D * M::RowsAtCompileTime;
    static constexpr int cols = M::ColsAtCompileTime == Eigen::Dynamic
                                    ? Eigen::Dynamic
                                    : D * M::ColsAtCompileTime;
    using Scalar = typename M::Scalar;
    using Matrix = Eigen::Matrix<Scalar, rows, cols, Eigen::ColMajor>;

    constexpr explicit EigenBlockMatrixViewFunctor(const M& matrix)
        : matrix_(matrix){};

    constexpr const Scalar operator()(Eigen::Index row, Eigen::Index col) const
    {
        if (row / matrix_.rows() != col / matrix_.cols())
        {
            return 0;
        }
        return matrix_(row % matrix_.rows(), col % matrix_.cols());
    }

private:
    const typename M::Nested& matrix_;
};

template <int D, typename M>
constexpr Eigen::CwiseNullaryOp<
    EigenBlockMatrixViewFunctor<D, M>,
    typename EigenBlockMatrixViewFunctor<D, M>::Matrix>
eigenBlockMatrixView(const Eigen::MatrixBase<M>& matrix)
{
    using Matrix = typename EigenBlockMatrixViewFunctor<D, M>::Matrix;
    return Matrix::NullaryExpr(
        D * matrix.rows(), D * matrix.cols(),
        EigenBlockMatrixViewFunctor<D, M>(matrix.derived()));
}
}  // namespace MathLib
