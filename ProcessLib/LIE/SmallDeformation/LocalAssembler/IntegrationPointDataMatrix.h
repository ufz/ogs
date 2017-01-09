/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESSLIB_LIE_SMALLDEFORMATION_INTEGRATIONPOINTDATAMATRIX_H_
#define PROCESSLIB_LIE_SMALLDEFORMATION_INTEGRATIONPOINTDATAMATRIX_H_

#include <memory>
#include <vector>

#include "MaterialLib/SolidModels/MechanicsBase.h"

namespace ProcessLib
{
namespace LIE
{
namespace SmallDeformation
{

template <typename BMatricesType, int DisplacementDim>
struct IntegrationPointDataMatrix final
{
    explicit IntegrationPointDataMatrix(
        MaterialLib::Solids::MechanicsBase<DisplacementDim>& solid_material)
        : _solid_material(solid_material),
          _material_state_variables(
              _solid_material.createMaterialStateVariables())
    {
    }

#if defined(_MSC_VER) && _MSC_VER < 1900
    // The default generated move-ctor is correctly generated for other
    // compilers.
    explicit IntegrationPointDataMatrix(IntegrationPointDataMatrix&& other)
        : _b_matrices(std::move(other._b_matrices)),
          _sigma(std::move(other._sigma)),
          _sigma_prev(std::move(other._sigma_prev)),
          _eps(std::move(other._eps)),
          _eps_prev(std::move(other._eps_prev)),
          _solid_material(other._solid_material),
          _material_state_variables(std::move(other._material_state_variables)),
          _C(std::move(other._C)),
          _detJ(std::move(other._detJ))
    {
    }
#endif  // _MSC_VER

    typename BMatricesType::BMatrixType _b_matrices;
    typename BMatricesType::KelvinVectorType _sigma, _sigma_prev;
    typename BMatricesType::KelvinVectorType _eps, _eps_prev;

    MaterialLib::Solids::MechanicsBase<DisplacementDim>& _solid_material;
    std::unique_ptr<
        typename MaterialLib::Solids::MechanicsBase<DisplacementDim>::MaterialStateVariables>
        _material_state_variables;

    typename BMatricesType::KelvinMatrixType _C;
    double _detJ;
    double _integralMeasure;

    void pushBackState()
    {
        _eps_prev = _eps;
        _sigma_prev = _sigma;
        _material_state_variables->pushBackState();
    }
};

}  // namespace SmallDeformation
}  // namespace LIE
}  // namespace ProcessLib

#endif // PROCESSLIB_LIE_SMALLDEFORMATION_INTEGRATIONPOINTDATAMATRIX_H_
