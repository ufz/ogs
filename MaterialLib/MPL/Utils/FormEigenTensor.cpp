/*
 * \file
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
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
        MaterialPropertyLib::Vector const& values) const
    {
        return Eigen::Map<Eigen::Matrix<double, GlobalDim, 1> const>(
                   values.data(), GlobalDim, 1)
            .asDiagonal();
    }

    Eigen::Matrix<double, GlobalDim, GlobalDim> operator()(
        MaterialPropertyLib::Tensor2d const& values) const
    {
        return Eigen::Map<Eigen::Matrix<double, GlobalDim, GlobalDim> const>(
            values.data(), GlobalDim, GlobalDim);
    }

    Eigen::Matrix<double, GlobalDim, GlobalDim> operator()(
        MaterialPropertyLib::Tensor const& values) const
    {
        return Eigen::Map<Eigen::Matrix<double, GlobalDim, GlobalDim> const>(
            values.data(), GlobalDim, GlobalDim);
    }

    Eigen::Matrix<double, GlobalDim, GlobalDim> operator()(
        MaterialPropertyLib::SymmTensor const& /*values*/) const
    {
        OGS_FATAL(
            "The value of MaterialPropertyLib::SymmTensor is inapplicable");
    }

    Eigen::Matrix<double, GlobalDim, GlobalDim> operator()(
        std::string const& /*values*/) const
    {
        OGS_FATAL("The value of std::string is inapplicable");
    }

    Eigen::Matrix<double, GlobalDim, GlobalDim> operator()(
        MaterialPropertyLib::Pair const& /*values*/) const
    {
        OGS_FATAL("The size of tensor is neither one nor %d nor %d squared.",
                  GlobalDim, GlobalDim);
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
