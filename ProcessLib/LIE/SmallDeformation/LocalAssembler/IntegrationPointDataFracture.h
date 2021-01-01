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

#include <Eigen/Eigen>

#include "MaterialLib/FractureModels/FractureModelBase.h"

namespace ProcessLib
{
namespace LIE
{
namespace SmallDeformation
{
template <typename HMatricesType, int DisplacementDim>
struct IntegrationPointDataFracture final
{
    explicit IntegrationPointDataFracture(
        MaterialLib::Fracture::FractureModelBase<DisplacementDim>&
            fracture_material_)
        : fracture_material(fracture_material_),
          material_state_variables(
              fracture_material.createMaterialStateVariables())
    {
    }

    typename HMatricesType::HMatrixType h_matrices;
    typename HMatricesType::ForceVectorType sigma, sigma_prev;
    typename HMatricesType::ForceVectorType w, w_prev;
    double aperture = 0.0;
    double aperture_prev = 0.0;
    double aperture0 = 0.0;

    MaterialLib::Fracture::FractureModelBase<DisplacementDim>&
        fracture_material;
    std::unique_ptr<typename MaterialLib::Fracture::FractureModelBase<
        DisplacementDim>::MaterialStateVariables>
        material_state_variables;

    Eigen::MatrixXd C;
    double integration_weight;

    void pushBackState()
    {
        w_prev = w;
        sigma_prev = sigma;
        aperture_prev = aperture;
        material_state_variables->pushBackState();
    }

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

}  // namespace SmallDeformation
}  // namespace LIE
}  // namespace ProcessLib
