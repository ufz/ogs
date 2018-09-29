/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
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
            fracture_material)
        : _fracture_material(fracture_material),
          _material_state_variables(
              _fracture_material.createMaterialStateVariables())
    {
    }

    typename HMatricesType::HMatrixType _h_matrices;
    typename HMatricesType::ForceVectorType _sigma, _sigma_prev;
    typename HMatricesType::ForceVectorType _w, _w_prev;
    double _aperture = 0.0;
    double _aperture_prev = 0.0;
    double _aperture0 = 0.0;

    MaterialLib::Fracture::FractureModelBase<DisplacementDim>&
        _fracture_material;
    std::unique_ptr<typename MaterialLib::Fracture::FractureModelBase<
        DisplacementDim>::MaterialStateVariables>
        _material_state_variables;

    Eigen::MatrixXd _C;
    double integration_weight;

    void pushBackState()
    {
        _w_prev = _w;
        _sigma_prev = _sigma;
        _aperture_prev = _aperture;
        _material_state_variables->pushBackState();
    }

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

}  // namespace SmallDeformation
}  // namespace LIE
}  // namespace ProcessLib
