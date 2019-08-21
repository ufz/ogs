/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  \file   StaggeredHTFEM-impl.h
 *  Created on October 13, 2017, 3:52 PM
 */

#pragma once

#include "StaggeredHTFEM.h"

#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/MPL/Utils/FormEigenTensor.h"
#include "ProcessLib/CoupledSolutionsForStaggeredScheme.h"

namespace ProcessLib
{
namespace HT
{
template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
void StaggeredHTFEM<ShapeFunction, IntegrationMethod, GlobalDim>::
    assembleForStaggeredScheme(double const t,
                               std::vector<double>& local_M_data,
                               std::vector<double>& local_K_data,
                               std::vector<double>& local_b_data,
                               LocalCoupledSolutions const& coupled_xs)
{
    if (coupled_xs.process_id == _heat_transport_process_id)
    {
        assembleHeatTransportEquation(t, local_M_data, local_K_data,
                                      local_b_data, coupled_xs);
        return;
    }

    assembleHydraulicEquation(t, local_M_data, local_K_data, local_b_data,
                              coupled_xs);
}

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
void StaggeredHTFEM<ShapeFunction, IntegrationMethod, GlobalDim>::
    assembleHydraulicEquation(double const t, std::vector<double>& local_M_data,
                              std::vector<double>& local_K_data,
                              std::vector<double>& local_b_data,
                              LocalCoupledSolutions const& coupled_xs)
{
    auto const& local_p = coupled_xs.local_coupled_xs[_hydraulic_process_id];
    auto const local_matrix_size = local_p.size();
    // This assertion is valid only if all nodal d.o.f. use the same shape
    // matrices.
    assert(local_matrix_size == ShapeFunction::NPOINTS);

    auto const& local_T1 =
        coupled_xs.local_coupled_xs[_heat_transport_process_id];
    auto const& local_T0 =
        coupled_xs.local_coupled_xs0[_heat_transport_process_id];
    const double dt = coupled_xs.dt;

    auto local_M = MathLib::createZeroedMatrix<LocalMatrixType>(
        local_M_data, local_matrix_size, local_matrix_size);
    auto local_K = MathLib::createZeroedMatrix<LocalMatrixType>(
        local_K_data, local_matrix_size, local_matrix_size);
    auto local_b = MathLib::createZeroedVector<LocalVectorType>(
        local_b_data, local_matrix_size);

    ParameterLib::SpatialPosition pos;
    pos.setElementID(this->_element.getID());

    auto const& process_data = this->_process_data;
    auto const& medium = *this->_process_data.media_map->getMedium(
        this->_element.getID());
    auto const& liquid_phase = medium.phase("AqueousLiquid");
    auto const& solid_phase = medium.phase("Solid");

    auto const& b = process_data.specific_body_force;

    MaterialPropertyLib::VariableArray vars;

    unsigned const n_integration_points =
        this->_integration_method.getNumberOfPoints();

    for (unsigned ip(0); ip < n_integration_points; ip++)
    {
        pos.setIntegrationPoint(ip);

        auto const& ip_data = this->_ip_data[ip];
        auto const& N = ip_data.N;
        auto const& dNdx = ip_data.dNdx;
        auto const& w = ip_data.integration_weight;

        double p_int_pt = 0.0;
        double T1_int_pt = 0.0;
        NumLib::shapeFunctionInterpolate(local_p, N, p_int_pt);
        NumLib::shapeFunctionInterpolate(local_T1, N, T1_int_pt);

        vars[static_cast<int>(MaterialPropertyLib::Variable::temperature)] =
            T1_int_pt;
        vars[static_cast<int>(MaterialPropertyLib::Variable::phase_pressure)] =
            p_int_pt;

        auto const porosity =
            solid_phase.property(MaterialPropertyLib::PropertyType::porosity)
                .template value<double>(vars, pos, t);
        auto const fluid_density =
            liquid_phase.property(MaterialPropertyLib::PropertyType::density)
                .template value<double>(vars, pos, t);

        const double dfluid_density_dp =
            liquid_phase.property(MaterialPropertyLib::PropertyType::density)
                .template dValue<double>(
                    vars, MaterialPropertyLib::Variable::phase_pressure, pos,
                    t);

        // Use the viscosity model to compute the viscosity
        auto const viscosity =
            liquid_phase.property(MaterialPropertyLib::PropertyType::viscosity)
                .template value<double>(vars, pos, t);

        // \todo the argument to getValue() has to be changed for non
        // constant storage model
        auto const specific_storage =
            solid_phase.property(MaterialPropertyLib::PropertyType::storage)
                .template value<double>(vars, pos, t);

        auto const intrinsic_permeability =
            MaterialPropertyLib::formEigenTensor<GlobalDim>(
                solid_phase
                    .property(MaterialPropertyLib::PropertyType::permeability)
                    .value(vars, pos, t));
        GlobalDimMatrixType const K_over_mu =
            intrinsic_permeability / viscosity;

        // matrix assembly
        local_M.noalias() +=
            w *
            (porosity * dfluid_density_dp / fluid_density + specific_storage) *
            N.transpose() * N;

        local_K.noalias() += w * dNdx.transpose() * K_over_mu * dNdx;

        if (process_data.has_gravity)
        {
            local_b.noalias() +=
                w * fluid_density * dNdx.transpose() * K_over_mu * b;
        }

        if (!process_data.has_fluid_thermal_expansion)
        {
            return;
        }

        // Add the thermal expansion term
        {
            auto const solid_thermal_expansion =
                process_data.solid_thermal_expansion(t, pos)[0];
            const double dfluid_density_dT =
                liquid_phase
                    .property(MaterialPropertyLib::PropertyType::density)
                    .template dValue<double>(
                        vars, MaterialPropertyLib::Variable::temperature, pos,
                        t);
            double T0_int_pt = 0.;
            NumLib::shapeFunctionInterpolate(local_T0, N, T0_int_pt);
            auto const biot_constant =
                process_data.biot_constant(t, pos)[0];
            const double eff_thermal_expansion =
                3.0 * (biot_constant - porosity) * solid_thermal_expansion -
                porosity * dfluid_density_dT / fluid_density;
            local_b.noalias() +=
                eff_thermal_expansion * (T1_int_pt - T0_int_pt) * w * N / dt;
        }
    }
}

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
void StaggeredHTFEM<ShapeFunction, IntegrationMethod, GlobalDim>::
    assembleHeatTransportEquation(double const t,
                                  std::vector<double>& local_M_data,
                                  std::vector<double>& local_K_data,
                                  std::vector<double>& /*local_b_data*/,
                                  LocalCoupledSolutions const& coupled_xs)
{
    auto const& local_p = coupled_xs.local_coupled_xs[_hydraulic_process_id];
    auto const local_matrix_size = local_p.size();
    // This assertion is valid only if all nodal d.o.f. use the same shape
    // matrices.
    assert(local_matrix_size == ShapeFunction::NPOINTS);

    auto local_p_Eigen_type =
        Eigen::Map<const NodalVectorType>(&local_p[0], local_matrix_size);

    auto const& local_T1 =
        coupled_xs.local_coupled_xs[_heat_transport_process_id];

    auto local_M = MathLib::createZeroedMatrix<LocalMatrixType>(
        local_M_data, local_matrix_size, local_matrix_size);
    auto local_K = MathLib::createZeroedMatrix<LocalMatrixType>(
        local_K_data, local_matrix_size, local_matrix_size);

    ParameterLib::SpatialPosition pos;
    pos.setElementID(this->_element.getID());

    auto const& process_data = this->_process_data;
    auto const& medium =
        *process_data.media_map->getMedium(this->_element.getID());
    auto const& liquid_phase = medium.phase("AqueousLiquid");
    auto const& solid_phase = medium.phase("Solid");

    auto const& b = process_data.specific_body_force;

    GlobalDimMatrixType const& I(
        GlobalDimMatrixType::Identity(GlobalDim, GlobalDim));

    MaterialPropertyLib::VariableArray vars;

    unsigned const n_integration_points =
        this->_integration_method.getNumberOfPoints();

    for (unsigned ip(0); ip < n_integration_points; ip++)
    {
        pos.setIntegrationPoint(ip);

        auto const& ip_data = this->_ip_data[ip];
        auto const& N = ip_data.N;
        auto const& dNdx = ip_data.dNdx;
        auto const& w = ip_data.integration_weight;

        double p_at_xi = 0.;
        NumLib::shapeFunctionInterpolate(local_p, N, p_at_xi);
        double T1_at_xi = 0.;
        NumLib::shapeFunctionInterpolate(local_T1, N, T1_at_xi);

        vars[static_cast<int>(MaterialPropertyLib::Variable::temperature)] =
            T1_at_xi;
        vars[static_cast<int>(MaterialPropertyLib::Variable::phase_pressure)] =
            p_at_xi;

        auto const porosity =
            solid_phase.property(MaterialPropertyLib::PropertyType::porosity)
                .template value<double>(vars, pos, t);

        // Use the fluid density model to compute the density
        auto const fluid_density =
            liquid_phase.property(MaterialPropertyLib::PropertyType::density)
                .template value<double>(vars, pos, t);
        auto const specific_heat_capacity_fluid =
            liquid_phase.property(MaterialPropertyLib::specific_heat_capacity)
                .template value<double>(vars, pos, t);

        // Assemble mass matrix
        local_M.noalias() += w *
                             this->getHeatEnergyCoefficient(
                                 vars, porosity, fluid_density,
                                 specific_heat_capacity_fluid, pos, t) *
                             N.transpose() * N;

        // Assemble Laplace matrix
        auto const viscosity =
            liquid_phase.property(MaterialPropertyLib::PropertyType::viscosity)
                .template value<double>(vars, pos, t);

        auto const intrinsic_permeability =
            MaterialPropertyLib::formEigenTensor<GlobalDim>(
                solid_phase
                    .property(MaterialPropertyLib::PropertyType::permeability)
                    .value(vars, pos, t));

        GlobalDimMatrixType const K_over_mu =
            intrinsic_permeability / viscosity;
        GlobalDimVectorType const velocity =
            process_data.has_gravity
                ? GlobalDimVectorType(-K_over_mu * (dNdx * local_p_Eigen_type -
                                                    fluid_density * b))
                : GlobalDimVectorType(-K_over_mu * dNdx * local_p_Eigen_type);

        GlobalDimMatrixType const thermal_conductivity_dispersivity =
            this->getThermalConductivityDispersivity(
                vars, porosity, fluid_density, specific_heat_capacity_fluid,
                velocity, I, pos, t);

        local_K.noalias() +=
            w * (dNdx.transpose() * thermal_conductivity_dispersivity * dNdx +
                 N.transpose() * velocity.transpose() * dNdx * fluid_density *
                     specific_heat_capacity_fluid);
    }
}

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
std::vector<double> const&
StaggeredHTFEM<ShapeFunction, IntegrationMethod, GlobalDim>::
    getIntPtDarcyVelocity(const double t,
                          GlobalVector const& /*current_solution*/,
                          NumLib::LocalToGlobalIndexMap const& dof_table,
                          std::vector<double>& cache) const
{
    auto const indices = NumLib::getIndices(this->_element.getID(), dof_table);
    assert(!indices.empty());
    std::vector<std::vector<GlobalIndexType>> indices_of_all_coupled_processes =
        {indices, indices};
    auto const local_xs = getCurrentLocalSolutions(
        *(this->_coupled_solutions), indices_of_all_coupled_processes);

    return this->getIntPtDarcyVelocityLocal(
        t, local_xs[_hydraulic_process_id],
        local_xs[_heat_transport_process_id], cache);
}
}  // namespace HT
}  // namespace ProcessLib
