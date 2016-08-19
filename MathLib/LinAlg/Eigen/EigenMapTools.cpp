/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "EigenMapTools.h"
#include <cassert>

namespace MathLib
{
Eigen::Map<
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>
toZeroedMatrix(std::vector<double>& data,
               Eigen::MatrixXd::Index rows,
               Eigen::MatrixXd::Index cols)
{
    assert(data.empty());  // in order that resize fills the vector with zeros.
    data.resize(rows * cols);
    return {data.data(), rows, cols};
}

Eigen::Map<const Eigen::
               Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>
toMatrix(std::vector<double> const& data,
         Eigen::MatrixXd::Index rows,
         Eigen::MatrixXd::Index cols)
{
    assert(static_cast<Eigen::MatrixXd::Index>(data.size()) == rows * cols);
    return {data.data(), rows, cols};
}

Eigen::Map<
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>
toMatrix(std::vector<double>& data,
         Eigen::MatrixXd::Index rows,
         Eigen::MatrixXd::Index cols)
{
    assert(static_cast<Eigen::MatrixXd::Index>(data.size()) == rows * cols);
    return {data.data(), rows, cols};
}

Eigen::Map<Eigen::VectorXd> toZeroedVector(std::vector<double>& data,
                                           Eigen::VectorXd::Index rows)
{
    assert(data.empty());  // in order that resize fills the vector with zeros.
    data.resize(rows);
    return {data.data(), rows};
}

Eigen::Map<const Eigen::VectorXd> toVector(std::vector<double> const& data,
                                           Eigen::VectorXd::Index rows)
{
    assert(static_cast<Eigen::VectorXd::Index>(data.size()) == rows);
    return {data.data(), rows};
}

Eigen::Map<Eigen::VectorXd> toVector(std::vector<double>& data,
                                     Eigen::VectorXd::Index rows)
{
    assert(static_cast<Eigen::VectorXd::Index>(data.size()) == rows);
    return {data.data(), rows};
}
}  // MathLib
