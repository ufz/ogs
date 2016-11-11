/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESSLIB_LIE_HYDROMECHANICS_INTEGRATIONPOINTDATAMATRIX_H_
#define PROCESSLIB_LIE_HYDROMECHANICS_INTEGRATIONPOINTDATAMATRIX_H_

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
          typename ShapeMatrixTypePressure, unsigned GlobalDim, unsigned NPoints>
struct IntegrationPointDataMatrix final
{
    explicit IntegrationPointDataMatrix(
        MaterialLib::Solids::MechanicsBase<GlobalDim>& solid_material)
        : solid_material(solid_material),
          material_state_variables(
              solid_material.createMaterialStateVariables())
    {
    }

    // The default generated move-ctor is correctly generated for other
    // compilers.
    explicit IntegrationPointDataMatrix(IntegrationPointDataMatrix&& other)
        : N_u(std::move(other.N_u)),
          H_u(std::move(other.H_u)),
          b_matrices(std::move(other.b_matrices)),
          sigma_eff(std::move(other.sigma_eff)),
          sigma_eff_prev(std::move(other.sigma_eff_prev)),
          eps(std::move(other.eps)),
          eps_prev(std::move(other.eps_prev)),
          N_p(std::move(other.N_p)),
          dNdx_p(std::move(other.dNdx_p)),
          solid_material(other.solid_material),
          material_state_variables(std::move(other.material_state_variables)),
          C(std::move(other.C)),
          integration_weight(std::move(other.integration_weight)),
          darcy_velocity(std::move(other.darcy_velocity))
    {
    }

    typename ShapeMatrixTypeDisplacement::NodalRowVectorType N_u;
    typename ShapeMatrixTypeDisplacement::template MatrixType<
        GlobalDim, NPoints * GlobalDim>
        H_u;
    typename BMatricesType::BMatrixType b_matrices;
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

    Eigen::Vector3d darcy_velocity;

    void pushBackState()
    {
        eps_prev = eps;
        sigma_eff_prev = sigma_eff;
        material_state_variables->pushBackState();
    }
};

}  // namespace HydroMechanics
}  // namespace LIE
}  // namespace ProcessLib

#endif // PROCESSLIB_LIE_HYDROMECHANICS_INTEGRATIONPOINTDATAMATRIX_H_
