/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
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
          typename ShapeMatrixTypePressure, int GlobalDim, unsigned NPoints>
struct IntegrationPointDataMatrix final
{
    explicit IntegrationPointDataMatrix(
        MaterialLib::Solids::MechanicsBase<GlobalDim>& solid_material)
        : solid_material(solid_material),
          material_state_variables(
              solid_material.createMaterialStateVariables()),
          darcy_velocity(GlobalDimVector::Zero())
    {
    }

    using GlobalDimVector = Eigen::Matrix<double, GlobalDim, 1>;

    typename ShapeMatrixTypeDisplacement::NodalRowVectorType N_u;
    typename ShapeMatrixTypeDisplacement::GlobalDimNodalMatrixType dNdx_u;
    typename ShapeMatrixTypeDisplacement::template MatrixType<
        GlobalDim, NPoints * GlobalDim>
        H_u;
    typename BMatricesType::KelvinVectorType sigma_eff, sigma_eff_prev;
    typename BMatricesType::KelvinVectorType eps, eps_prev;

    typename ShapeMatrixTypePressure::NodalRowVectorType N_p;
    typename ShapeMatrixTypePressure::GlobalDimNodalMatrixType dNdx_p;

    MaterialLib::Solids::MechanicsBase<GlobalDim>& solid_material;
    std::unique_ptr<typename MaterialLib::Solids::MechanicsBase<
        GlobalDim>::MaterialStateVariables>
        material_state_variables;

    typename BMatricesType::KelvinMatrixType C;
    double integration_weight;

    GlobalDimVector darcy_velocity;

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
