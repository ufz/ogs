/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
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
          typename ShapeMatrixTypePressure, int GlobalDim>
struct IntegrationPointDataFracture final
{
    explicit IntegrationPointDataFracture(
        MaterialLib::Fracture::FractureModelBase<GlobalDim>& fracture_material_)
        : fracture_material(fracture_material_),
          material_state_variables(
              fracture_material.createMaterialStateVariables()),
          darcy_velocity(GlobalDimVectorType::Zero())
    {
    }

    using GlobalDimVectorType = Eigen::Matrix<double, GlobalDim, 1>;

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

    std::unique_ptr<
        typename MaterialLib::Fracture::Permeability::PermeabilityState>
        permeability_state;

    Eigen::MatrixXd C;
    double integration_weight;

    GlobalDimVectorType darcy_velocity;

    void pushBackState()
    {
        w_prev = w;
        sigma_eff_prev = sigma_eff;
        material_state_variables->pushBackState();
    }

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

}  // namespace HydroMechanics
}  // namespace LIE
}  // namespace ProcessLib
