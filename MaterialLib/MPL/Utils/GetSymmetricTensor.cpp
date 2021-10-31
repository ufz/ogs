/*
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on March 04, 2020, 5:20 PM
 */

#include "GetSymmetricTensor.h"

#include "MaterialLib/MPL/PropertyType.h"

namespace MaterialPropertyLib
{
template <int GlobalDim>
struct GetSymmetricTensor
{
    SymmetricTensor<GlobalDim> operator()(double const& value) const
    {
        SymmetricTensor<GlobalDim> result = SymmetricTensor<GlobalDim>::Zero();
        result.template head<3>() = Eigen::Vector3d::Constant(value);
        return result;
    }

    SymmetricTensor<GlobalDim> operator()(Eigen::Vector2d const& values) const
    {
        SymmetricTensor<GlobalDim> result = SymmetricTensor<GlobalDim>::Zero();
        result.template head<2>() = values;
        return result;
    }

    SymmetricTensor<GlobalDim> operator()(Eigen::Vector3d const& values) const
    {
        SymmetricTensor<GlobalDim> result = SymmetricTensor<GlobalDim>::Zero();
        result.template head<3>() = values;
        return result;
    }

    SymmetricTensor<GlobalDim> operator()(Eigen::Matrix2d const& values) const
    {
        if constexpr (GlobalDim == 2)
        {
            SymmetricTensor<GlobalDim> result;
            result << values(0, 0), values(1, 1), 0., values(0, 1);
            return result;
        }
        OGS_FATAL("Cannot convert 2d matrix to 3d symmetric Tensor.");
    }

    SymmetricTensor<GlobalDim> operator()(Eigen::Matrix3d const& values) const
    {
        if constexpr (GlobalDim == 3)
        {
            SymmetricTensor<GlobalDim> result;
            result << values(0, 0), values(1, 1), values(2, 2), values(0, 1),
                values(1, 2), values(0, 2);
            return result;
        }
        OGS_FATAL("Cannot convert 3d matrix to 2d symmetric Tensor.");
    }

    SymmetricTensor<GlobalDim> operator()(
        SymmetricTensor<2> const& values) const
    {
        if constexpr (GlobalDim == 2)
        {
            return values;
        }
        OGS_FATAL("Cannot convert 3d symmetric tensor to 2d symmetric tensor.");
    }

    SymmetricTensor<GlobalDim> operator()(
        SymmetricTensor<3> const& values) const
    {
        if constexpr (GlobalDim == 3)
        {
            return values;
        }
        OGS_FATAL("Cannot convert 2d symmetric tensor to 3d symmetric tensor.");
    }

    SymmetricTensor<GlobalDim> operator()(
        Eigen::MatrixXd const& /*values*/) const
    {
        OGS_FATAL(
            "Cannot convert dynamic Eigen matrix to {:d}d symmetric tensor.",
            GlobalDim);
    }
};

template <int GlobalDim>
SymmetricTensor<GlobalDim> getSymmetricTensor(
    MaterialPropertyLib::PropertyDataType const& values)
{
    return std::visit(GetSymmetricTensor<GlobalDim>(), values);
}

template Eigen::Matrix<double, 4, 1> getSymmetricTensor<2>(
    MaterialPropertyLib::PropertyDataType const& values);

template Eigen::Matrix<double, 6, 1> getSymmetricTensor<3>(
    MaterialPropertyLib::PropertyDataType const& values);

}  // namespace MaterialPropertyLib
