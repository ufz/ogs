/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "MaterialLib/MPL/Properties/EmbeddedFracturePermeability.h"

#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/MPL/Utils/FormEigenTensor.h"

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
    double const threshold_strain,
    ParameterLib::Parameter<double> const& fracture_rotation_xy,
    ParameterLib::Parameter<double> const& fracture_rotation_yz,
    double const jacobian_factor)
    : _n(fracture_normal),
      _n_const(fracture_normal_is_constant),
      _k(intrinsic_permeability),
      _b0(initial_aperture),
      _a(mean_fracture_distance),
      _e0(threshold_strain),
      _phi_xy(fracture_rotation_xy),
      _phi_yz(fracture_rotation_yz),
      _jf(jacobian_factor)
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
    ParameterLib::SpatialPosition const& pos, double const t,
    double const /*dt*/) const
{
    Eigen::Matrix<double, 3, 1> const n = [&]
    {
        if (_n_const)
        {
            return _n;
        }
        auto const sigma = formEigenTensor<3>(std::get<SymmetricTensor>(
            variable_array[static_cast<int>(Variable::total_stress)]));
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, 3, 3>> e_s(sigma);
        return (Eigen::Matrix<double, 3, 1>)e_s.eigenvectors().col(2);
    }();

    auto const rotMat_xy =
        Eigen::AngleAxisd(_phi_xy(t, pos)[0], Eigen::Vector3d::UnitZ());
    auto const rotMat_yz =
        Eigen::AngleAxisd(_phi_yz(t, pos)[0], Eigen::Vector3d::UnitX());

    Eigen::Matrix<double, 3, 1> const n_r = rotMat_yz * (rotMat_xy * n);

    auto const eps = MathLib::KelvinVector::kelvinVectorToTensor(
        std::get<MathLib::KelvinVector::KelvinVectorType<DisplacementDim>>(
            variable_array[static_cast<int>(Variable::mechanical_strain)]));
    double const e_n = (eps * n_r).dot(n_r.transpose());
    double const H_de = (e_n > _e0) ? 1.0 : 0.0;
    double const b_f = _b0 + H_de * _a * (e_n - _e0);
    double const coeff = H_de * (b_f / _a) * ((b_f * b_f / 12.0) - _k);

    Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
    return (coeff * (I - n_r * n_r.transpose()) + _k * I)
        .template topLeftCorner<DisplacementDim, DisplacementDim>()
        .eval();
}

template <int DisplacementDim>
PropertyDataType EmbeddedFracturePermeability<DisplacementDim>::dValue(
    VariableArray const& variable_array, Variable const primary_variable,
    ParameterLib::SpatialPosition const& pos, double const t,
    double const /*dt*/) const
{
    if (primary_variable != Variable::mechanical_strain)
    {
        OGS_FATAL(
            "EmbeddedFracturePermeability::dValue is implemented for "
            "derivatives with respect to strain only.");
    }

    Eigen::Matrix<double, 3, 1> const n = [&]
    {
        if (_n_const)
        {
            return _n;
        }
        auto const sigma = formEigenTensor<3>(std::get<SymmetricTensor>(
            variable_array[static_cast<int>(Variable::total_stress)]));
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, 3, 3>> e_s(sigma);
        return (Eigen::Matrix<double, 3, 1>)e_s.eigenvectors().col(2);
    }();

    auto const rotMat_xy =
        Eigen::AngleAxisd(_phi_xy(t, pos)[0], Eigen::Vector3d::UnitZ());
    auto const rotMat_yz =
        Eigen::AngleAxisd(_phi_yz(t, pos)[0], Eigen::Vector3d::UnitX());

    Eigen::Matrix<double, 3, 1> const n_r = rotMat_yz * (rotMat_xy * n);

    auto const eps = MathLib::KelvinVector::kelvinVectorToTensor(
        std::get<MathLib::KelvinVector::KelvinVectorType<DisplacementDim>>(
            variable_array[static_cast<int>(Variable::mechanical_strain)]));
    double const e_n = (eps * n_r).dot(n_r.transpose());
    double const H_de = (e_n > _e0) ? 1.0 : 0.0;
    double const b_f = _b0 + H_de * _a * (e_n - _e0);

    Eigen::Matrix3d const M = n_r * n_r.transpose();
    return Eigen::MatrixXd(
        _jf * H_de * (b_f * b_f / 4 - _k) *
        MathLib::KelvinVector::tensorToKelvin<DisplacementDim>(
            Eigen::Matrix3d::Identity() - M) *
        MathLib::KelvinVector::tensorToKelvin<DisplacementDim>(M).transpose());
}

template class EmbeddedFracturePermeability<2>;
template class EmbeddedFracturePermeability<3>;
}  // namespace MaterialPropertyLib
