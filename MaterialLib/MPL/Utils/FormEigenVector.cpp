/*
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on July 31, 2019, 11:28 AM
 */

#include "FormEigenVector.h"

#include "MaterialLib/MPL/PropertyType.h"

namespace MaterialPropertyLib
{
template <int GlobalDim>
struct FormEigenVector
{
    Eigen::Matrix<double, GlobalDim, 1> operator()(double const value) const
    {
        if constexpr (GlobalDim == 1)
        {
            return Eigen::Matrix<double, 1, 1>{value};
        }
        if constexpr (GlobalDim == 2)
        {
            return Eigen::Matrix<double, 2, 1>{value, value};
        }
        if constexpr (GlobalDim == 3)
        {
            return Eigen::Matrix<double, 3, 1>{value, value, value};
        }
    }

    Eigen::Matrix<double, GlobalDim, 1> operator()(
        Eigen::Vector2d const& values) const
    {
        if constexpr (GlobalDim == 2)
        {
            return values;
        }
        OGS_FATAL("Cannot convert 2d vector to {:d}d vector.", GlobalDim);
    }

    Eigen::Matrix<double, GlobalDim, 1> operator()(
        Eigen::Vector3d const& values) const
    {
        if constexpr (GlobalDim == 3)
        {
            return values;
        }
        OGS_FATAL("Cannot convert 3d vector to a {:d}d vector.", GlobalDim);
    }

    Eigen::Matrix<double, GlobalDim, 1> operator()(
        Eigen::Matrix<double, 2, 2> const& /*values*/) const
    {
        OGS_FATAL("Cannot convert a 2d tensor to a {:d}d Vector.", GlobalDim);
    }
    Eigen::Matrix<double, GlobalDim, 1> operator()(
        Eigen::Matrix<double, 3, 3> const& /*values*/) const
    {
        OGS_FATAL("Cannot convert a 3d tensor to a {:d}d Vector.", GlobalDim);
    }

    Eigen::Matrix<double, GlobalDim, 1> operator()(
        Eigen::Matrix<double, 4, 1> const& /*values*/) const
    {
        OGS_FATAL("Cannot convert a 4d vector to a {:d}d vector.", GlobalDim);
    }

    Eigen::Matrix<double, GlobalDim, 1> operator()(
        Eigen::Matrix<double, 6, 1> const& /*values*/) const
    {
        OGS_FATAL("Cannot convert a 6d vector to a {:d}d vector.", GlobalDim);
    }

    Eigen::Matrix<double, GlobalDim, 1> operator()(
        Eigen::MatrixXd const& /*values*/) const
    {
        OGS_FATAL("Cannot convert dynamic Eigen matrix to a {:d}d vector ",
                  GlobalDim);
    }
};

template <int GlobalDim>
Eigen::Matrix<double, GlobalDim, 1> formEigenVector(
    MaterialPropertyLib::PropertyDataType const& values)
{
    return std::visit(FormEigenVector<GlobalDim>(), values);
}

template Eigen::Matrix<double, 1, 1> formEigenVector<1>(
    MaterialPropertyLib::PropertyDataType const& values);

template Eigen::Matrix<double, 2, 1> formEigenVector<2>(
    MaterialPropertyLib::PropertyDataType const& values);

template Eigen::Matrix<double, 3, 1> formEigenVector<3>(
    MaterialPropertyLib::PropertyDataType const& values);

}  // namespace MaterialPropertyLib
