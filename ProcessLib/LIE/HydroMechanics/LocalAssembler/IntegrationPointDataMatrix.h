/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <memory>
#include <vector>

#include "MaterialLib/SolidModels/MechanicsBase.h"

namespace ProcessLib
{
namespace LIE
{
namespace HydroMechanics
{
template <typename BMatricesType, typename ShapeMatrixTypeDisplacement,
          typename ShapeMatrixTypePressure, int DisplacementDim,
          unsigned NPoints>
struct IntegrationPointDataMatrix final
{
    explicit IntegrationPointDataMatrix(
        MaterialLib::Solids::MechanicsBase<DisplacementDim>& solid_material)
        : solid_material(solid_material),
          material_state_variables(
              solid_material.createMaterialStateVariables()),
          darcy_velocity(DisplacementDimVector::Zero())
    {
    }

    using DisplacementDimVector = Eigen::Matrix<double, DisplacementDim, 1>;

    typename ShapeMatrixTypeDisplacement::NodalRowVectorType N_u;
    typename ShapeMatrixTypeDisplacement::GlobalDimNodalMatrixType dNdx_u;
    typename ShapeMatrixTypeDisplacement::template MatrixType<
        DisplacementDim, NPoints * DisplacementDim>
        H_u;
    typename BMatricesType::KelvinVectorType sigma_eff, sigma_eff_prev;
    typename BMatricesType::KelvinVectorType eps, eps_prev;

    typename ShapeMatrixTypePressure::NodalRowVectorType N_p;
    typename ShapeMatrixTypePressure::GlobalDimNodalMatrixType dNdx_p;

    using GlobalDimVectorType =
        typename ShapeMatrixTypePressure::GlobalDimVectorType;

    MaterialLib::Solids::MechanicsBase<DisplacementDim>& solid_material;
    std::unique_ptr<typename MaterialLib::Solids::MechanicsBase<
        DisplacementDim>::MaterialStateVariables>
        material_state_variables;

    typename BMatricesType::KelvinMatrixType C;
    double integration_weight;

    GlobalDimVectorType darcy_velocity;

    void pushBackState()
    {
        eps_prev = eps;
        sigma_eff_prev = sigma_eff;
        material_state_variables->pushBackState();
    }

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

}  // namespace HydroMechanics
}  // namespace LIE
}  // namespace ProcessLib
