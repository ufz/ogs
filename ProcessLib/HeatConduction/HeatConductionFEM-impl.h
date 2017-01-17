/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  \file   HeatConductionFEM-impl.h
 *  Created on January 17, 2017, 3:41 PM
 */

#pragma once

#include "HeatConductionFEM.h"

namespace ProcessLib
{
namespace HeatConduction
{
template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
template <typename LiquidFlowVelocityCalculator>
void LocalAssemblerData<ShapeFunction, IntegrationMethod, GlobalDim>::
    assembleHeatTransportLiquidFlow(
        double const t, int const material_id, SpatialPosition& pos,
        int const gravitational_axis_id,
        double const gravitational_acceleration, Eigen::MatrixXd const& perm,
        ProcessLib::LiquidFlow::LiquidFlowMaterialProperties const&
            liquid_flow_prop,
        std::vector<double> const& local_x, std::vector<double> const& local_p,
        std::vector<double>& local_M_data, std::vector<double>& local_K_data)
{
    auto const local_matrix_size = local_x.size();
    // This assertion is valid only if all nodal d.o.f. use the same shape
    // matrices.
    assert(local_matrix_size == ShapeFunction::NPOINTS * NUM_NODAL_DOF);

    auto local_M = MathLib::createZeroedMatrix<NodalMatrixType>(
        local_M_data, local_matrix_size, local_matrix_size);
    auto local_K = MathLib::createZeroedMatrix<NodalMatrixType>(
        local_K_data, local_matrix_size, local_matrix_size);
    const auto local_p_vec =
        MathLib::toVector<NodalVectorType>(local_p, local_matrix_size);

    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    const double porosity_variable = 0.;
    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        pos.setIntegrationPoint(ip);
        auto const& sm = _shape_matrices[ip];
        auto const& wp = _integration_method.getWeightedPoint(ip);
        double p = 0.;
        NumLib::shapeFunctionInterpolate(local_p, sm.N, p);
        double T = 0.;
        NumLib::shapeFunctionInterpolate(local_x, sm.N, T);

        // Material parameters of solid phase
        auto const k_s = _process_data.thermal_conductivity(t, pos)[0];
        auto const cp_s = _process_data.heat_capacity(t, pos)[0];
        auto const rho_s = _process_data.density(t, pos)[0];

        // Material parameters of fluid phase
        double const cp_f = liquid_flow_prop.getHeatCapacity(p, T);
        double const k_f = liquid_flow_prop.getThermalConductivity(p, T);
        double const rho_f = liquid_flow_prop.getLiquidDensity(p, T);

        // Material parameter of porosity
        double const poro =
            liquid_flow_prop.getPorosity(material_id, porosity_variable, T);

        double const effective_cp =
            (1.0 - poro) * cp_s * rho_s + poro * cp_f * rho_f;
        double const effective_K = (1.0 - poro) * k_s + poro * k_f;

        double const integration_factor =
            sm.detJ * wp.getWeight() * sm.integralMeasure;

        local_M.noalias() +=
            effective_cp * sm.N.transpose() * sm.N * integration_factor;
        local_K.noalias() +=
            sm.dNdx.transpose() * effective_K * sm.dNdx * integration_factor;

        // Compute fluid flow velocity
        double const mu = liquid_flow_prop.getViscosity(p, T);  // Viscosity
        GlobalDimVectorType const& darcy_velocity =
            LiquidFlowVelocityCalculator::calculateVelocity(
                local_p_vec, sm, perm, mu, rho_f * gravitational_acceleration,
                gravitational_axis_id);
        // Add the advection term
        local_K.noalias() += cp_f * rho_f * sm.N.transpose() *
                             (darcy_velocity.transpose() * sm.dNdx) *
                             integration_factor;
    }
}

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
void LocalAssemblerData<ShapeFunction, IntegrationMethod, GlobalDim>::
    coupling_assemble(double const t, std::vector<double> const& local_x,
                      std::vector<double>& local_M_data,
                      std::vector<double>& local_K_data,
                      std::vector<double>& /*local_b_data*/,
                      LocalCouplingTerm const& coupled_term)
{
    auto it = coupled_term.coupled_processes.begin();
    while (it != coupled_term.coupled_processes.end())
    {
        switch (it->first)
        {
            case ProcessLib::ProcessType::LiquidFlowProcess:
            {
                ProcessLib::LiquidFlow::LiquidFlowProcess const& pcs =
                    static_cast<
                        ProcessLib::LiquidFlow::LiquidFlowProcess const&>(
                        it->second);
                const auto liquid_flow_prop =
                    pcs.getLiquidFlowMaterialProperties();

                const auto local_p = coupled_term.local_coupled_xs.at(
                    ProcessLib::ProcessType::LiquidFlowProcess);

                SpatialPosition pos;
                pos.setElementID(_element.getID());
                const int material_id = liquid_flow_prop->getMaterialID(pos);

                const Eigen::MatrixXd& perm = liquid_flow_prop->getPermeability(
                    material_id, t, pos, _element.getDimension());

                if (perm.size() == 1)
                {
                    assembleHeatTransportLiquidFlow<
                        IsotropicLiquidFlowVelocityCalculator>(
                        t, material_id, pos, pcs.getGravitationalAxisID(),
                        pcs.getGravitationalacceleration(), perm,
                        *liquid_flow_prop, local_x, local_p, local_M_data,
                        local_K_data);
                }
                else
                {
                    assembleHeatTransportLiquidFlow<
                        AnisotropicLiquidFlowVelocityCalculator>(
                        t, material_id, pos, pcs.getGravitationalAxisID(),
                        pcs.getGravitationalacceleration(), perm,
                        *liquid_flow_prop, local_x, local_p, local_M_data,
                        local_K_data);
                }
            }
            break;
            default:
                OGS_FATAL(
                    "This coupled process is not presented for "
                    "HeatConduction process");
        }
        it++;
    }
}

}  // namespace HeatConduction
}  // namespace ProcessLib
