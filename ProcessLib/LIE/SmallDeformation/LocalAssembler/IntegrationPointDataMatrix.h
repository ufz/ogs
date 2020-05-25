/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
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
namespace SmallDeformation
{
template <typename ShapeMatricesType, typename BMatricesType,
          int DisplacementDim>
struct IntegrationPointDataMatrix final
{
    explicit IntegrationPointDataMatrix(
        MaterialLib::Solids::MechanicsBase<DisplacementDim>& solid_material)
        : solid_material_(solid_material),
          material_state_variables_(
              solid_material_.createMaterialStateVariables())
    {
    }

    typename BMatricesType::KelvinVectorType sigma_, sigma_prev_;
    typename BMatricesType::KelvinVectorType eps_, eps_prev_;

    MaterialLib::Solids::MechanicsBase<DisplacementDim>& solid_material_;
    std::unique_ptr<typename MaterialLib::Solids::MechanicsBase<
        DisplacementDim>::MaterialStateVariables>
        material_state_variables_;

    typename BMatricesType::KelvinMatrixType C_;
    double integration_weight;

    typename ShapeMatricesType::NodalRowVectorType N;
    typename ShapeMatricesType::GlobalDimNodalMatrixType dNdx;

    void pushBackState()
    {
        eps_prev_ = eps_;
        sigma_prev_ = sigma_;
        material_state_variables_->pushBackState();
    }

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

}  // namespace SmallDeformation
}  // namespace LIE
}  // namespace ProcessLib
