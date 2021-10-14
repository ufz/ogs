/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "MaterialLib/MPL/Properties/OrthotropicEmbeddedFracturePermeability.h"

#include "MaterialLib/MPL/Medium.h"

namespace MaterialPropertyLib
{
template <int DisplacementDim>
OrthotropicEmbeddedFracturePermeability<DisplacementDim>::
    OrthotropicEmbeddedFracturePermeability(
        std::string name,
        std::vector<double> const& mean_fracture_distances,
        std::vector<double> const& threshold_strains,
        Eigen::Matrix<double, 3, 3> const fracture_normals,
        ParameterLib::Parameter<double> const& intrinsic_permeability,
        ParameterLib::Parameter<double> const& fracture_rotation_xy,
        ParameterLib::Parameter<double> const& fracture_rotation_yz,
        double const jacobian_factor)
    : _a(mean_fracture_distances),
      _e0(threshold_strains),
      _n(fracture_normals),
      _k(intrinsic_permeability),
      _phi_xy(fracture_rotation_xy),
      _phi_yz(fracture_rotation_yz),
      _jf(jacobian_factor)
{
    name_ = std::move(name);
}

template <int DisplacementDim>
void OrthotropicEmbeddedFracturePermeability<DisplacementDim>::checkScale()
    const
{
    if (!std::holds_alternative<Medium*>(scale_))
    {
        OGS_FATAL(
            "The property 'OrthotropicEmbeddedFracturePermeability' is "
            "implemented on the 'media' scale only.");
    }
}

template <int DisplacementDim>
PropertyDataType
OrthotropicEmbeddedFracturePermeability<DisplacementDim>::value(
    VariableArray const& variable_array,
    ParameterLib::SpatialPosition const& pos, double const t,
    double const /*dt*/) const
{
    auto const eps = MathLib::KelvinVector::kelvinVectorToTensor(
        std::get<MathLib::KelvinVector::KelvinVectorType<DisplacementDim>>(
            variable_array[static_cast<int>(Variable::mechanical_strain)]));
    double const k = std::get<double>(fromVector(_k(t, pos)));
    double const _b0 = std::sqrt(12.0 * k);

    double const phi_xy = std::get<double>(fromVector(_phi_xy(t, pos)));
    double const phi_yz = std::get<double>(fromVector(_phi_yz(t, pos)));
    auto const rotMat_xy = Eigen::AngleAxisd(phi_xy, Eigen::Vector3d::UnitZ());
    auto const rotMat_yz = Eigen::AngleAxisd(phi_yz, Eigen::Vector3d::UnitX());

    Eigen::Matrix3d result = Eigen::Vector3d::Constant(k).asDiagonal();

    for (int i = 0; i < 3; i++)
    {
        Eigen::Vector3d const n_i = rotMat_yz * (rotMat_xy * _n.col(i));
        double const e_n = (eps * n_i).dot(n_i.transpose());
        double const H_de = (e_n > _e0[i]) ? 1.0 : 0.0;
        double const b_f = _b0 + H_de * _a[i] * (e_n - _e0[i]);

        // The H_de factor is only valid as long as _b0 equals
        // std::sqrt(12.0 * k), else it has to be dropped.
        result += H_de * (b_f / _a[i]) * ((b_f * b_f / 12.0) - k) *
                  (Eigen::Matrix3d::Identity() - n_i * n_i.transpose());
    }

    return result.template topLeftCorner<DisplacementDim, DisplacementDim>()
        .eval();
}
template <int DisplacementDim>
PropertyDataType
OrthotropicEmbeddedFracturePermeability<DisplacementDim>::dValue(
    VariableArray const& variable_array, Variable const primary_variable,
    ParameterLib::SpatialPosition const& pos, double const t,
    double const /*dt*/) const
{
    if (primary_variable != Variable::mechanical_strain)
    {
        OGS_FATAL(
            "OrthotropicEmbeddedFracturePermeability::dValue is implemented "
            "for derivatives with respect to strain only.");
    }

    auto const eps = MathLib::KelvinVector::kelvinVectorToTensor(
        std::get<MathLib::KelvinVector::KelvinVectorType<DisplacementDim>>(
            variable_array[static_cast<int>(Variable::mechanical_strain)]));
    double const k = std::get<double>(fromVector(_k(t, pos)));
    double const _b0 = std::sqrt(12.0 * k);

    double const phi_xy = std::get<double>(fromVector(_phi_xy(t, pos)));
    double const phi_yz = std::get<double>(fromVector(_phi_yz(t, pos)));
    auto const rotMat_xy = Eigen::AngleAxisd(phi_xy, Eigen::Vector3d::UnitZ());
    auto const rotMat_yz = Eigen::AngleAxisd(phi_yz, Eigen::Vector3d::UnitX());

    MathLib::KelvinVector::KelvinMatrixType<DisplacementDim> result =
        MathLib::KelvinVector::KelvinMatrixType<DisplacementDim>::Zero();

    for (int i = 0; i < 3; i++)
    {
        Eigen::Vector3d const n_i = rotMat_yz * (rotMat_xy * _n.col(i));
        Eigen::Matrix3d const M = n_i * n_i.transpose();
        double const e_n = (eps * n_i).dot(n_i.transpose());
        double const H_de = (e_n > _e0[i]) ? 1.0 : 0.0;
        double const b_f = _b0 + H_de * _a[i] * (e_n - _e0[i]);

        result += H_de * (b_f * b_f / 4.0 - k) *
                  MathLib::KelvinVector::tensorToKelvin<DisplacementDim>(
                      Eigen::Matrix3d::Identity() - M) *
                  MathLib::KelvinVector::tensorToKelvin<DisplacementDim>(M)
                      .transpose();
    }
    return Eigen::MatrixXd(_jf * result);
}

template class OrthotropicEmbeddedFracturePermeability<2>;
template class OrthotropicEmbeddedFracturePermeability<3>;
}  // namespace MaterialPropertyLib
