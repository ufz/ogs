/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <Eigen/Dense>
#include <array>
#include <vector>

namespace ParameterLib
{
template <typename T>
struct Parameter;
class SpatialPosition;
}  // namespace ParameterLib

namespace ParameterLib
{
struct CoordinateSystem final
{
    CoordinateSystem(Parameter<double> const& e0, Parameter<double> const& e1);

    CoordinateSystem(Parameter<double> const& e0,
                     Parameter<double> const& e1,
                     Parameter<double> const& e2);

    template <int Dimension>
    Eigen::Matrix<double, Dimension, Dimension> transformation(
        SpatialPosition const& pos) const;

    Eigen::Matrix<double, 3, 3> transformation_3d(
        SpatialPosition const& pos) const;

    template <int Dimension>
    Eigen::Matrix<double, Dimension, Dimension> rotateTensor(
        std::vector<double> const& values, SpatialPosition const& pos) const;

    template <int Dimension>
    Eigen::Matrix<double, Dimension, Dimension> rotateDiagonalTensor(
        std::vector<double> const& values, SpatialPosition const& pos) const;

private:
    std::array<Parameter<double> const*, 3> _base;
};

extern template Eigen::Matrix<double, 2, 2> CoordinateSystem::rotateTensor<2>(
    std::vector<double> const& values, SpatialPosition const& pos) const;
extern template Eigen::Matrix<double, 3, 3> CoordinateSystem::rotateTensor<3>(
    std::vector<double> const& values, SpatialPosition const& pos) const;
extern template Eigen::Matrix<double, 2, 2>
CoordinateSystem::rotateDiagonalTensor<2>(std::vector<double> const& values,
                                          SpatialPosition const& pos) const;
extern template Eigen::Matrix<double, 3, 3>
CoordinateSystem::rotateDiagonalTensor<3>(std::vector<double> const& values,
                                          SpatialPosition const& pos) const;
}  // namespace ParameterLib
