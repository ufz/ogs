/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <vector>
#include <Eigen/Core>

namespace MathLib
{
Eigen::Map<
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>
toZeroedMatrix(std::vector<double>& data,
               Eigen::MatrixXd::Index rows,
               Eigen::MatrixXd::Index cols);

Eigen::Map<const Eigen::
               Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>
toMatrix(std::vector<double> const& data,
         Eigen::MatrixXd::Index rows,
         Eigen::MatrixXd::Index cols);

Eigen::Map<
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>
toMatrix(std::vector<double>& data,
         Eigen::MatrixXd::Index rows,
         Eigen::MatrixXd::Index cols);

Eigen::Map<Eigen::VectorXd> toZeroedVector(std::vector<double>& data,
                                           Eigen::VectorXd::Index rows);

Eigen::Map<const Eigen::VectorXd> toVector(std::vector<double> const& data,
                                           Eigen::VectorXd::Index rows);

Eigen::Map<Eigen::VectorXd> toVector(std::vector<double>& data,
                                     Eigen::VectorXd::Index rows);
}  // MathLib
