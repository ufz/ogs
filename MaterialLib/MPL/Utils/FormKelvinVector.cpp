/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on February 12, 2021, 4:34 PM
 */

#include "FormKelvinVector.h"

#include "BaseLib/Error.h"

namespace MaterialPropertyLib
{
static constexpr const char error_info[] =
    "The conversion to a Kelvin vector of correct dimensionality is ambiguous."
    "Please use a scalar number for isotropic properties, a three element "
    "array or a 3 x 3 matrix for anisotropic properties.";

template <int GlobalDim>
struct FormKelvinVector
{
    MathLib::KelvinVector::KelvinVectorType<GlobalDim> operator()(
        double const& value) const
    {
        MathLib::KelvinVector::KelvinVectorType<GlobalDim> result =
            MathLib::KelvinVector::KelvinVectorType<GlobalDim>::Zero();
        result.template head<3>() = Eigen::Vector3d::Constant(value);
        return result;
    }

    MathLib::KelvinVector::KelvinVectorType<GlobalDim> operator()(
        Eigen::Matrix<double, 2, 1> const& /*values*/) const
    {
        OGS_FATAL(error_info);
    }

    MathLib::KelvinVector::KelvinVectorType<GlobalDim> operator()(
        Eigen::Matrix<double, 3, 1> const& values) const
    {
        MathLib::KelvinVector::KelvinVectorType<GlobalDim> result =
            MathLib::KelvinVector::KelvinVectorType<GlobalDim>::Zero();
        result.template head<3>() = values;
        return result;
    }

    MathLib::KelvinVector::KelvinVectorType<GlobalDim> operator()(
        Eigen::Matrix<double, 2, 2> const& /*values*/) const
    {
        OGS_FATAL(error_info);
    }

    MathLib::KelvinVector::KelvinVectorType<GlobalDim> operator()(
        Eigen::Matrix<double, 3, 3> const& values) const
    {
        return MathLib::KelvinVector::tensorToKelvin<GlobalDim>(values);
    }

    MathLib::KelvinVector::KelvinVectorType<GlobalDim> operator()(
        Eigen::Matrix<double, 4, 1> const& /*values*/) const
    {
        OGS_FATAL(error_info);
    }

    MathLib::KelvinVector::KelvinVectorType<GlobalDim> operator()(
        Eigen::Matrix<double, 6, 1> const& /*values*/) const
    {
        OGS_FATAL(error_info);
    }

    MathLib::KelvinVector::KelvinVectorType<GlobalDim> operator()(
        Eigen::MatrixXd const& /*values*/) const
    {
        OGS_FATAL(error_info);
    }
};

template <int GlobalDim>
MathLib::KelvinVector::KelvinVectorType<GlobalDim> formKelvinVector(
    MaterialPropertyLib::PropertyDataType const& values)
{
    return std::visit(FormKelvinVector<GlobalDim>(), values);
}

template MathLib::KelvinVector::KelvinVectorType<2> formKelvinVector<2>(
    MaterialPropertyLib::PropertyDataType const& values);

template MathLib::KelvinVector::KelvinVectorType<3> formKelvinVector<3>(
    MaterialPropertyLib::PropertyDataType const& values);

}  // namespace MaterialPropertyLib
