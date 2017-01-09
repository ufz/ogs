/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <Eigen/Eigen>

#include "MaterialLib/FractureModels/FractureModelBase.h"

namespace ProcessLib
{
namespace LIE
{
namespace HydroMechanics
{
template <typename HMatricesType, typename ShapeMatrixTypeDisplacement,
          typename ShapeMatrixTypePressure, unsigned GlobalDim>
struct IntegrationPointDataFracture final
{
    explicit IntegrationPointDataFracture(
        MaterialLib::Fracture::FractureModelBase<GlobalDim>& fracture_material_)
        : fracture_material(fracture_material_),
          material_state_variables(
              fracture_material.createMaterialStateVariables())
    {
    }

    IntegrationPointDataFracture(IntegrationPointDataFracture&& other)
        : H_u(std::move(other.H_u)),
          sigma_eff(std::move(other.sigma_eff)),
          sigma_eff_prev(std::move(other.sigma_eff_prev)),
          w(std::move(other.w)),
          w_prev(std::move(other.w_prev)),
          N_p(std::move(other.N_p)),
          dNdx_p(std::move(other.dNdx_p)),
          aperture(std::move(other.aperture)),
          aperture0(std::move(other.aperture0)),
          fracture_material(other.fracture_material),
          C(std::move(other.C)),
          integration_weight(std::move(other.integration_weight)),
          darcy_velocity(std::move(other.darcy_velocity))
    {
    }

    typename HMatricesType::HMatrixType H_u;
    typename HMatricesType::ForceVectorType sigma_eff, sigma_eff_prev;
    typename HMatricesType::ForceVectorType w, w_prev;

    typename ShapeMatrixTypePressure::NodalRowVectorType N_p;
    typename ShapeMatrixTypePressure::GlobalDimNodalMatrixType dNdx_p;

    double aperture = 0.0;
    double aperture0 = 0.0;
    double permeability = 0.0;

    MaterialLib::Fracture::FractureModelBase<GlobalDim>& fracture_material;
    std::unique_ptr<typename MaterialLib::Fracture::FractureModelBase<
        GlobalDim>::MaterialStateVariables>
        material_state_variables;

    Eigen::MatrixXd C;
    double integration_weight;

    Eigen::Vector3d darcy_velocity;

    void pushBackState()
    {
        w_prev = w;
        sigma_eff_prev = sigma_eff;
    }
};

}  // namespace HydroMechanics
}  // namespace LIE
}  // namespace ProcessLib
