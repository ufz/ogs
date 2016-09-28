/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file   LiquidFlowLocalAssembler.cpp
 *
 * Created on August 19, 2016, 2:28 PM
 */

#include "LiquidFlowLocalAssembler.h"

#include "NumLib/Function/Interpolation.h"
#include "MaterialLib/PhysicalConstant.h"

namespace ProcessLib
{
namespace LiquidFlow
{
template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
void LiquidFlowLocalAssembler<ShapeFunction, IntegrationMethod, GlobalDim>::
    assemble(double const /*t*/, std::vector<double> const& local_p,
             std::vector<double>& local_M_data,
             std::vector<double>& local_K_data,
             std::vector<double>& local_b_data)
{
    auto const local_matrix_size = local_p.size();
    assert(local_matrix_size == ShapeFunction::NPOINTS);

    auto local_M = MathLib::createZeroedMatrix<NodalMatrixType>(
        local_M_data, local_matrix_size, local_matrix_size);
    auto local_K = MathLib::createZeroedMatrix<NodalMatrixType>(
        local_K_data, local_matrix_size, local_matrix_size);
    auto local_b = MathLib::createZeroedVector<NodalVectorType>(
        local_b_data, local_matrix_size);

    // For velocity calculation:
    // TODO: move velocity calculation to a member function, and call the
    // function after the linearization step if the calculation is needed.
    const auto local_p_vec =
        MathLib::toVector<NodalVectorType>(local_p, local_matrix_size);

    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    const unsigned mat_id = 0; // TODO for heterogeneous medium
    const MaterialLib::PorousMedium::CoefMatrix& perm =
        _material_properties.intrinsic_permeabiliy[mat_id];

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        auto const& sm = _shape_matrices[ip];
        auto const& wp = _integration_method.getWeightedPoint(ip);

        double p = 0.;
        NumLib::shapeFunctionInterpolate(local_p, sm.N, p);
        // TODO : compute _temperature from the heat transport pcs

        // Assemble mass matrix, M
        local_M.noalias() +=
            _material_properties.getMassCoeffcient(p, _temperature, mat_id) *
            sm.N.transpose() * sm.N * sm.detJ * wp.getWeight();

        // Compute density:
        const double rho_g =
            _material_properties.getLiquidDensity(p, _temperature) *
            MaterialLib::PhysicalConstant::g;
        // Compute viscosity:
        const double mu = _material_properties.getViscosity(p, _temperature);

        // Assemble Laplacian, K, and RHS by the gravitational term
        if (perm.size() == 1)  // Save time for isotropic permeability.
        {
            //  Use scalar number for isotropic permeability
            //  to save the computation time.
            const double K = perm(0, 0) / mu;
            const double fac = K * sm.detJ * wp.getWeight();
            local_K.noalias() += fac * sm.dNdx.transpose() * sm.dNdx;
            local_b.noalias() -=
                fac * sm.dNdx.transpose().col(GlobalDim - 1) * rho_g;

            // Compute the velocity
            GlobalDimVectorType darcy_velocity = -K * sm.dNdx * local_p_vec;
            // gravity term
            darcy_velocity[GlobalDim - 1] -= K * rho_g;
            for (unsigned d = 0; d < GlobalDim; ++d)
            {
                _darcy_velocities[d][ip] = darcy_velocity[d];
            }
        }
        else
        {
            const double fac = sm.detJ * wp.getWeight() / mu;
            local_K.noalias() += fac * sm.dNdx.transpose() * perm * sm.dNdx;
            local_b.noalias() -=
                fac * rho_g * sm.dNdx.transpose() * perm.col(GlobalDim - 1);

            // Compute the velocity
            GlobalDimVectorType const darcy_velocity =
                -perm * sm.dNdx * local_p_vec / mu -
                rho_g * perm.col(GlobalDim - 1) / mu;
            for (unsigned d = 0; d < GlobalDim; ++d)
            {
                _darcy_velocities[d][ip] = darcy_velocity[d];
            }
        }
    }
}

}  // end of namespace
}  // end of namespace
