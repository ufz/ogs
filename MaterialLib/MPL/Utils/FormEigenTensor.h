// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <Eigen/Core>
#include <type_traits>
#include <variant>

#include "MaterialLib/MPL/Property.h"
#include "MathLib/FormattingUtils.h"

namespace MaterialPropertyLib
{
template <int GlobalDim>
struct FormEigenTensor
{
    constexpr Eigen::Matrix<double, GlobalDim, GlobalDim> operator()(
        double const& value) const
    {
        return Eigen::Matrix<double, GlobalDim, GlobalDim>::Identity() * value;
    }

    // Constexpr version for GlobalDim == 2 (valid case)
    template <int Dim = GlobalDim>
    constexpr std::enable_if_t<Dim == 2, Eigen::Matrix<double, 2, 2>>
    operator()(Eigen::Vector2d const& values) const
    {
        return values.asDiagonal();
    }

    // Non-constexpr version for invalid GlobalDim values
    template <int Dim = GlobalDim>
    constexpr std::enable_if_t<Dim != 2,
                               Eigen::Matrix<double, GlobalDim, GlobalDim>>
    operator()(Eigen::Vector2d const& values) const
    {
        OGS_FATAL(
            "Cannot convert 2d vector with values [{}] to {:d}x{:d} diagonal "
            "matrix.",
            values, GlobalDim, GlobalDim);
    }

    // Constexpr version for GlobalDim == 3 (valid case)
    template <int Dim = GlobalDim>
    constexpr std::enable_if_t<Dim == 3, Eigen::Matrix<double, 3, 3>>
    operator()(Eigen::Vector3d const& values) const
    {
        return values.asDiagonal();
    }

    // Non-constexpr version for invalid GlobalDim values
    template <int Dim = GlobalDim>
    constexpr std::enable_if_t<Dim != 3,
                               Eigen::Matrix<double, GlobalDim, GlobalDim>>
    operator()(Eigen::Vector3d const& values) const
    {
        OGS_FATAL(
            "Cannot convert 3d vector with values [{}] to {:d}x{:d} diagonal "
            "matrix.",
            values, GlobalDim, GlobalDim);
    }

    // Constexpr version for GlobalDim == 2 (valid case)
    template <int Dim = GlobalDim>
    constexpr std::enable_if_t<Dim == 2, Eigen::Matrix<double, 2, 2>>
    operator()(Eigen::Matrix<double, 2, 2> const& values) const
    {
        return values;
    }

    // Non-constexpr version for invalid GlobalDim values
    template <int Dim = GlobalDim>
    constexpr std::enable_if_t<Dim != 2,
                               Eigen::Matrix<double, GlobalDim, GlobalDim>>
    operator()(Eigen::Matrix<double, 2, 2> const& values) const
    {
        OGS_FATAL(
            "Cannot convert a 2d tensor with values [{}] to {:d}x{:d} matrix",
            values, GlobalDim, GlobalDim);
    }

    // Constexpr version for GlobalDim == 3 (valid case)
    template <int Dim = GlobalDim>
    constexpr std::enable_if_t<Dim == 3, Eigen::Matrix<double, 3, 3>>
    operator()(Eigen::Matrix<double, 3, 3> const& values) const
    {
        return values;
    }

    // Non-constexpr version for invalid GlobalDim values
    template <int Dim = GlobalDim>
    constexpr std::enable_if_t<Dim != 3,
                               Eigen::Matrix<double, GlobalDim, GlobalDim>>
    operator()(Eigen::Matrix<double, 3, 3> const& values) const
    {
        OGS_FATAL(
            "Cannot convert a 3d tensor with values [{}] to {:d}x{:d} matrix",
            values, GlobalDim, GlobalDim);
    }

    constexpr Eigen::Matrix<double, GlobalDim, GlobalDim> operator()(
        Eigen::Matrix<double, 4, 1> const& values) const
    {
        if constexpr (GlobalDim == 2)
        {  // skip the z-direction in this case
            Eigen::Matrix<double, GlobalDim, GlobalDim> result;
            result << values[0], values[3], values[3], values[1];
            return result;
        }
        else if constexpr (GlobalDim == 3)
        {
            Eigen::Matrix<double, GlobalDim, GlobalDim> result;
            result << values[0], values[3], 0, values[3], values[1], 0, 0, 0,
                values[2];
            return result;
        }
        // For other dimensions, we need to provide a fallback to avoid missing
        // return warning
        OGS_FATAL("Cannot convert 4d vector to {:d}x{:d} matrix", GlobalDim,
                  GlobalDim);
    }

    // Constexpr version for GlobalDim == 3 (valid case)
    template <int Dim = GlobalDim>
    constexpr std::enable_if_t<Dim == 3, Eigen::Matrix<double, 3, 3>>
    operator()(Eigen::Matrix<double, 6, 1> const& values) const
    {
        Eigen::Matrix<double, 3, 3> result;
        result << values[0], values[3], values[5], values[3], values[1],
            values[4], values[5], values[4], values[2];
        return result;
    }

    // Non-constexpr version for invalid GlobalDim values
    template <int Dim = GlobalDim>
    constexpr std::enable_if_t<Dim != 3,
                               Eigen::Matrix<double, GlobalDim, GlobalDim>>
    operator()(Eigen::Matrix<double, 6, 1> const& values) const
    {
        OGS_FATAL(
            "Cannot convert a symmetric 3d tensor with values [{}] to a {}x{} "
            "matrix",
            values, GlobalDim, GlobalDim);
    }

    Eigen::Matrix<double, GlobalDim, GlobalDim> operator()(
        Eigen::MatrixXd const& values) const
    {
        if (GlobalDim == values.rows() && GlobalDim == values.cols())
        {
            return values;
        }

        OGS_FATAL(
            "Cannot convert a dynamic {}x{} matrix with values [{}] to a {}x{} "
            "matrix",
            values.rows(), values.cols(), values, GlobalDim, GlobalDim);
    }
};

template <int GlobalDim>
constexpr Eigen::Matrix<double, GlobalDim, GlobalDim> formEigenTensor(
    MaterialPropertyLib::PropertyDataType const& values)
{
    return std::visit(FormEigenTensor<GlobalDim>(), values);
}

}  // namespace MaterialPropertyLib
