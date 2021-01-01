/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include "IntegrationPointDataNonlocalInterface.h"
#include "MaterialLib/SolidModels/Ehlers.h"

namespace ProcessLib
{
namespace SmallDeformationNonlocal
{
template <typename BMatricesType, typename ShapeMatricesType,
          int DisplacementDim>
struct IntegrationPointData final : public IntegrationPointDataNonlocalInterface
{
    explicit IntegrationPointData(
        MaterialLib::Solids::Ehlers::SolidEhlers<DisplacementDim>&
            solid_material)
        : solid_material(solid_material),
          material_state_variables(
              solid_material.createMaterialStateVariables())
    {
        auto const msv =
            static_cast<typename MaterialLib::Solids::Ehlers::StateVariables<
                DisplacementDim>&>(*material_state_variables);

        eps_p_V = &msv.eps_p.V;
        eps_p_D_xx = &msv.eps_p.D[0];
    }

    typename BMatricesType::BMatrixType b_matrices;
    typename BMatricesType::KelvinVectorType sigma, sigma_prev;
    typename BMatricesType::KelvinVectorType eps, eps_prev;
    double free_energy_density = 0;
    double damage = 0;       ///< isotropic damage
    double damage_prev = 0;  ///< \copydoc damage
    double kappa_d_prev = 0;  ///< \copydoc kappa_d

    MaterialLib::Solids::Ehlers::SolidEhlers<DisplacementDim>& solid_material;
    std::unique_ptr<typename MaterialLib::Solids::MechanicsBase<
        DisplacementDim>::MaterialStateVariables>
        material_state_variables;

    typename BMatricesType::KelvinMatrixType C;
    typename ShapeMatricesType::NodalRowVectorType N;
    typename ShapeMatricesType::GlobalDimNodalMatrixType dNdx;

    double const* eps_p_V;     // Used for the secondary variables output.
    double const* eps_p_D_xx;  // Used for the secondary variables output.

    void pushBackState()
    {
        eps_prev = eps;
        sigma_prev = sigma;
        damage_prev = damage;
        kappa_d_prev = kappa_d;
        material_state_variables->pushBackState();
    }

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

}  // namespace SmallDeformationNonlocal
}  // namespace ProcessLib
