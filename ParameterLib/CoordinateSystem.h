/**
 * \file
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
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

    template <int Dimension>
    Eigen::Matrix<double, Dimension, Dimension> rotateTensor(
        std::vector<double> const& values, SpatialPosition const& pos) const
    {
        assert(values.size() == Dimension * Dimension ||
               "Input vector has wrong dimension; expected 4 or 9 entries.");
        auto const tensor =
            Eigen::Map<Eigen::Matrix<double, Dimension, Dimension> const>(
                values.data(), Dimension, Dimension);
        auto const R = transformation<Dimension>(pos);
        return R * tensor * R.transpose();
    }

    template <int Dimension>
    Eigen::Matrix<double, Dimension, Dimension> rotateDiagonalTensor(
        std::vector<double> const& values, SpatialPosition const& pos) const
    {
        assert(values.size() == Dimension ||
               "Input vector has wrong dimension; expected 2 or 3 entries.");
        auto const tensor =
            Eigen::Map<Eigen::Matrix<double, Dimension, 1> const>(values.data(),
                                                                  Dimension, 1);
        auto const R = transformation<Dimension>(pos);
        return R * tensor.asDiagonal() * R.transpose();
    }

private:
    std::array<Parameter<double> const*, 3> _base;
};

}  // namespace ParameterLib
