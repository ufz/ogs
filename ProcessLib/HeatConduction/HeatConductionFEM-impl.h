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
#include "NumLib/Function/Interpolation.h"

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
        double const gravitational_acceleration,
        Eigen::MatrixXd const& permeability,
        ProcessLib::LiquidFlow::LiquidFlowMaterialProperties const&
            liquid_flow_properties,
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

        // Thermal conductivity of solid phase.
        auto const k_s = _process_data.thermal_conductivity(t, pos)[0];
        // Specific heat capacity of solid phase.
        auto const cp_s = _process_data.heat_capacity(t, pos)[0];
        // Density of solid phase.
        auto const rho_s = _process_data.density(t, pos)[0];

        // Thermal conductivity of liquid.
        double const k_f = liquid_flow_properties.getThermalConductivity(p, T);
        // Specific heat capacity of liquid.
        double const cp_f = liquid_flow_properties.getHeatCapacity(p, T);
        // Density of liquid.
        double const rho_f = liquid_flow_properties.getLiquidDensity(p, T);

        // Porosity of porous media.
        double const n = liquid_flow_properties.getPorosity(
            material_id, porosity_variable, T);

        // Effective specific heat capacity.
        double const effective_cp = (1.0 - n) * cp_s * rho_s + n * cp_f * rho_f;
        // Effective thermal conductivity.
        double const effective_K = (1.0 - n) * k_s + n * k_f;

        double const integration_factor =
            sm.detJ * wp.getWeight() * sm.integralMeasure;

        local_M.noalias() +=
            effective_cp * sm.N.transpose() * sm.N * integration_factor;
        local_K.noalias() +=
            sm.dNdx.transpose() * effective_K * sm.dNdx * integration_factor;

        // Compute fluid flow velocity:
        double const mu =
            liquid_flow_properties.getViscosity(p, T);  // Viscosity
        GlobalDimVectorType const& darcy_velocity =
            LiquidFlowVelocityCalculator::calculateVelocity(
                local_p_vec, sm, permeability, mu,
                rho_f * gravitational_acceleration, gravitational_axis_id);
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
    for (auto const& coupled_process_pair : coupled_term.coupled_processes)
    {
        if (coupled_process_pair.first ==
            std::type_index(typeid(ProcessLib::LiquidFlow::LiquidFlowProcess)))
        {
            assert(
                dynamic_cast<const ProcessLib::LiquidFlow::LiquidFlowProcess*>(
                    &(coupled_process_pair.second)) != nullptr);

            ProcessLib::LiquidFlow::LiquidFlowProcess const& pcs =
                static_cast<ProcessLib::LiquidFlow::LiquidFlowProcess const&>(
                    coupled_process_pair.second);
            const auto liquid_flow_prop = pcs.getLiquidFlowMaterialProperties();

            const auto local_p =
                coupled_term.local_coupled_xs.at(coupled_process_pair.first);

            SpatialPosition pos;
            pos.setElementID(_element.getID());
            const int material_id = liquid_flow_prop->getMaterialID(pos);

            const Eigen::MatrixXd& K = liquid_flow_prop->getPermeability(
                material_id, t, pos, _element.getDimension());

            if (K.size() == 1)
            {
                assembleHeatTransportLiquidFlow<
                    IsotropicLiquidFlowVelocityCalculator>(
                    t, material_id, pos, pcs.getGravitationalAxisID(),
                    pcs.getGravitationalAcceleration(), K, *liquid_flow_prop,
                    local_x, local_p, local_M_data, local_K_data);
            }
            else
            {
                assembleHeatTransportLiquidFlow<
                    AnisotropicLiquidFlowVelocityCalculator>(
                    t, material_id, pos, pcs.getGravitationalAxisID(),
                    pcs.getGravitationalAcceleration(), K, *liquid_flow_prop,
                    local_x, local_p, local_M_data, local_K_data);
            }
        }
        else
        {
            OGS_FATAL(
                "This coupled process is not presented for "
                "HeatConduction process");
        }
    }
}

}  // namespace HeatConduction
}  // namespace ProcessLib
