/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "MaterialLib/MPL/Properties/EmbeddedFracturePermeability.h"

#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/MPL/Utils/FormEigenTensor.h"
#include "MathLib/KelvinVector.h"

namespace MaterialPropertyLib
{
template <int DisplacementDim>
EmbeddedFracturePermeability<DisplacementDim>::EmbeddedFracturePermeability(
    std::string name,
    Eigen::Matrix<double, 3, 1> const fracture_normal,
    bool const fracture_normal_is_constant,
    double const intrinsic_permeability,
    double const initial_aperture,
    double const mean_fracture_distance,
    double const threshold_strain)
    : _n(fracture_normal),
      _n_const(fracture_normal_is_constant),
      _k(intrinsic_permeability),
      _b0(initial_aperture),
      _a(mean_fracture_distance),
      _e0(threshold_strain)
{
    name_ = std::move(name);
}

template <int DisplacementDim>
void EmbeddedFracturePermeability<DisplacementDim>::checkScale() const
{
    if (!std::holds_alternative<Medium*>(scale_))
    {
        OGS_FATAL(
            "The property 'EmbeddedFracturePermeability' is "
            "implemented on the 'medium' scale only.");
    }
}

template <int DisplacementDim>
PropertyDataType EmbeddedFracturePermeability<DisplacementDim>::value(
    VariableArray const& variable_array,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    Eigen::Matrix<double, 3, 1> const n = [&] {
        if (_n_const)
        {
            return _n;
        }
        auto const sigma = formEigenTensor<3>(std::get<SymmetricTensor>(
            variable_array[static_cast<int>(Variable::total_stress)]));
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, 3, 3>> e_s(sigma);
        return (Eigen::Matrix<double, 3, 1>)e_s.eigenvectors().col(2);
    }();

    auto const eps = MathLib::KelvinVector::kelvinVectorToTensor(
        std::get<MathLib::KelvinVector::KelvinVectorType<DisplacementDim>>(
            variable_array[static_cast<int>(Variable::mechanical_strain)]));
    double const e_n = (eps * n).dot(n.transpose());
    double const H_de = (e_n > _e0) ? 1.0 : 0.0;
    double const b_f = _b0 + H_de * _a * (e_n - _e0);
    double const coeff = H_de * (b_f / _a) * ((b_f * b_f / 12.0) - _k);

    Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
    return (coeff * (I - n * n.transpose()) + _k * I).eval();
}

template <int DisplacementDim>
PropertyDataType EmbeddedFracturePermeability<DisplacementDim>::dValue(
    VariableArray const& variable_array, Variable const primary_variable,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    (void)primary_variable;
    assert((primary_variable == Variable::mechanical_strain) &&
           "EmbeddedFracturePermeability::dValue is implemented for "
           "derivatives with respect to strain only.");

    Eigen::Matrix<double, 3, 1> const n = [&] {
        if (_n_const)
        {
            return _n;
        }
        auto const sigma = formEigenTensor<3>(std::get<SymmetricTensor>(
            variable_array[static_cast<int>(Variable::total_stress)]));
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, 3, 3>> e_s(sigma);
        return (Eigen::Matrix<double, 3, 1>)e_s.eigenvectors().col(2);
    }();

    auto const eps = MathLib::KelvinVector::kelvinVectorToTensor(
        std::get<MathLib::KelvinVector::KelvinVectorType<DisplacementDim>>(
            variable_array[static_cast<int>(Variable::mechanical_strain)]));
    double const e_n = (eps * n).dot(n.transpose());
    double const H_de = (e_n > _e0) ? 1.0 : 0.0;
    double const b_f = _b0 + H_de * _a * (e_n - _e0);

    Eigen::Matrix3d const M = n * n.transpose();
    return (H_de * (b_f * b_f / 4 - _k) * (Eigen::Matrix3d::Identity() - M) * M)
        .eval();
}

template class EmbeddedFracturePermeability<2>;
template class EmbeddedFracturePermeability<3>;
}  // namespace MaterialPropertyLib
