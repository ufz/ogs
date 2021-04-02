/*
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on July 31, 2019, 11:28 AM
 */

#include "FormEigenTensor.h"

#include "MaterialLib/MPL/PropertyType.h"

namespace MaterialPropertyLib
{
template <int GlobalDim>
struct FormEigenTensor
{
    Eigen::Matrix<double, GlobalDim, GlobalDim> operator()(
        double const& value) const
    {
        return Eigen::Matrix<double, GlobalDim, GlobalDim>::Identity() * value;
    }

    Eigen::Matrix<double, GlobalDim, GlobalDim> operator()(
        Eigen::Vector2d const& values) const
    {
        if constexpr (GlobalDim == 2)
        {
            return values.asDiagonal();
        }
        OGS_FATAL("Cannot convert 2d vector to {:d}x{:d} diagonal matrix.",
                  GlobalDim, GlobalDim);
    }

    Eigen::Matrix<double, GlobalDim, GlobalDim> operator()(
        Eigen::Vector3d const& values) const
    {
        if constexpr (GlobalDim == 3)
        {
            return values.asDiagonal();
        }
        OGS_FATAL("Cannot convert 3d vector to {:d}x{:d} diagonal matrix.",
                  GlobalDim, GlobalDim);
    }

    Eigen::Matrix<double, GlobalDim, GlobalDim> operator()(
        Eigen::Matrix<double, 2, 2> const& values) const
    {
        if constexpr (GlobalDim == 2)
        {
            return values;
        }
        OGS_FATAL("Cannot convert a 2d tensor to {:d}x{:d} matrix", GlobalDim,
                  GlobalDim);
    }
    Eigen::Matrix<double, GlobalDim, GlobalDim> operator()(
        Eigen::Matrix<double, 3, 3> const& values) const
    {
        if constexpr (GlobalDim == 3)
        {
            return values;
        }
        OGS_FATAL("Cannot convert a 3d tensor to {:d}x{:d} matrix", GlobalDim,
                  GlobalDim);
    }

    Eigen::Matrix<double, GlobalDim, GlobalDim> operator()(
        Eigen::Matrix<double, 4, 1> const& values) const
    {
        Eigen::Matrix<double, GlobalDim, GlobalDim> result;
        if constexpr (GlobalDim == 2)
        {  // skip the z-direction in this case
            result << values[0], values[3], values[3], values[1];
        }
        if constexpr (GlobalDim == 3)
        {
            result << values[0], values[3], 0, values[3], values[1], 0, 0, 0,
                values[2];
        }
        return result;
    }

    Eigen::Matrix<double, GlobalDim, GlobalDim> operator()(
        Eigen::Matrix<double, 6, 1> const& values) const
    {
        if constexpr (GlobalDim == 3)
        {
            Eigen::Matrix<double, GlobalDim, GlobalDim> result;
            result << values[0], values[3], values[5], values[3], values[1],
                values[4], values[5], values[4], values[2];
            return result;
        }

        OGS_FATAL("Cannot convert a symmetric 3d tensor to {:d}x{:d} matrix",
                  GlobalDim);
    }
};

template <int GlobalDim>
Eigen::Matrix<double, GlobalDim, GlobalDim> formEigenTensor(
    MaterialPropertyLib::PropertyDataType const& values)
{
    return std::visit(FormEigenTensor<GlobalDim>(), values);
}

template Eigen::Matrix<double, 1, 1> formEigenTensor<1>(
    MaterialPropertyLib::PropertyDataType const& values);

template Eigen::Matrix<double, 2, 2> formEigenTensor<2>(
    MaterialPropertyLib::PropertyDataType const& values);

template Eigen::Matrix<double, 3, 3> formEigenTensor<3>(
    MaterialPropertyLib::PropertyDataType const& values);

}  // namespace MaterialPropertyLib
