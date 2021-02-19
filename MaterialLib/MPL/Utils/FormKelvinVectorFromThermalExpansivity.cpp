/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on February 12, 2021, 4:34 PM
 */

#include "FormKelvinVectorFromThermalExpansivity.h"

#include "BaseLib/Error.h"

namespace MaterialPropertyLib
{
static const char error_info[] =
    "The thermal expansivity can only be either a scalar number for "
    "isotropic thermal expansion or a three element array for anisotropic "
    "thermal expansion.";

template <int GlobalDim>
struct FormKelvinVectorFromThermalExpansivity
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
        Eigen::Matrix<double, 3, 3> const& /*values*/) const
    {
        OGS_FATAL(error_info);
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
};

template <int GlobalDim>
MathLib::KelvinVector::KelvinVectorType<GlobalDim>
formKelvinVectorFromThermalExpansivity(
    MaterialPropertyLib::PropertyDataType const& values)
{
    return std::visit(FormKelvinVectorFromThermalExpansivity<GlobalDim>(),
                      values);
}

template MathLib::KelvinVector::KelvinVectorType<2>
formKelvinVectorFromThermalExpansivity<2>(
    MaterialPropertyLib::PropertyDataType const& values);

template MathLib::KelvinVector::KelvinVectorType<3>
formKelvinVectorFromThermalExpansivity<3>(
    MaterialPropertyLib::PropertyDataType const& values);

}  // namespace MaterialPropertyLib
