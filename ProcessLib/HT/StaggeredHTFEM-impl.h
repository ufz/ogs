/**
 * \file
 * \copyright
 * Copyright (c) 2012-2023, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  Created on October 13, 2017, 3:52 PM
 */

#pragma once

#include <typeinfo>

#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/MPL/Utils/FormEigenTensor.h"
#include "ProcessLib/CoupledSolutionsForStaggeredScheme.h"
#include "StaggeredHTFEM.h"

namespace ProcessLib
{
namespace HT
{
template <typename ShapeFunction, int GlobalDim>
void StaggeredHTFEM<ShapeFunction, GlobalDim>::assembleForStaggeredScheme(
    double const t, double const dt, Eigen::VectorXd const& local_x,
    Eigen::VectorXd const& local_x_prev, int const process_id,
    std::vector<double>& local_M_data, std::vector<double>& local_K_data,
    std::vector<double>& local_b_data)
{
    if (process_id == this->_process_data.heat_transport_process_id)
    {
        assembleHeatTransportEquation(t, dt, local_x, local_M_data,
                                      local_K_data);
        return;
    }

    assembleHydraulicEquation(t, dt, local_x, local_x_prev, local_M_data,
                              local_K_data, local_b_data);
}

template <typename ShapeFunction, int GlobalDim>
void StaggeredHTFEM<ShapeFunction, GlobalDim>::assembleHydraulicEquation(
    double const t, double const dt, Eigen::VectorXd const& local_x,
    Eigen::VectorXd const& local_x_prev, std::vector<double>& local_M_data,
    std::vector<double>& local_K_data, std::vector<double>& local_b_data)
{
    auto const local_p =
        local_x.template segment<pressure_size>(pressure_index);
    auto const local_T =
        local_x.template segment<temperature_size>(temperature_index);

    auto const local_T_prev =
        local_x_prev.template segment<temperature_size>(temperature_index);

    auto local_M = MathLib::createZeroedMatrix<LocalMatrixType>(
        local_M_data, pressure_size, pressure_size);
    auto local_K = MathLib::createZeroedMatrix<LocalMatrixType>(
        local_K_data, pressure_size, pressure_size);
    auto local_b = MathLib::createZeroedVector<LocalVectorType>(local_b_data,
                                                                pressure_size);

    ParameterLib::SpatialPosition pos;
    pos.setElementID(this->_element.getID());

    auto const& process_data = this->_process_data;
    auto const& medium =
        *this->_process_data.media_map.getMedium(this->_element.getID());
    auto const& liquid_phase = medium.phase("AqueousLiquid");
    auto const& solid_phase = medium.phase("Solid");

    auto const& b =
        process_data
            .projected_specific_body_force_vectors[this->_element.getID()];

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
        double T_int_pt = 0.0;
        NumLib::shapeFunctionInterpolate(local_p, N, p_int_pt);
        NumLib::shapeFunctionInterpolate(local_T, N, T_int_pt);

        vars.temperature = T_int_pt;
        vars.phase_pressure = p_int_pt;

        vars.liquid_saturation = 1.0;

        auto const porosity =
            medium.property(MaterialPropertyLib::PropertyType::porosity)
                .template value<double>(vars, pos, t, dt);
        auto const fluid_density =
            liquid_phase.property(MaterialPropertyLib::PropertyType::density)
                .template value<double>(vars, pos, t, dt);

        vars.density = fluid_density;
        const double dfluid_density_dp =
            liquid_phase.property(MaterialPropertyLib::PropertyType::density)
                .template dValue<double>(
                    vars, MaterialPropertyLib::Variable::phase_pressure, pos, t,
                    dt);

        // Use the viscosity model to compute the viscosity
        auto const viscosity =
            liquid_phase.property(MaterialPropertyLib::PropertyType::viscosity)
                .template value<double>(vars, pos, t, dt);

        // \todo the argument to getValue() has to be changed for non
        // constant storage model
        auto const specific_storage =
            solid_phase.property(MaterialPropertyLib::PropertyType::storage)
                .template value<double>(vars, pos, t, dt);

        auto const intrinsic_permeability =
            MaterialPropertyLib::formEigenTensor<GlobalDim>(
                medium.property(MaterialPropertyLib::PropertyType::permeability)
                    .value(vars, pos, t, dt));
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
                        t, dt);
            double const Tdot_int_pt = (T_int_pt - local_T_prev.dot(N)) / dt;
            auto const biot_constant = process_data.biot_constant(t, pos)[0];
            const double eff_thermal_expansion =
                3.0 * (biot_constant - porosity) * solid_thermal_expansion -
                porosity * dfluid_density_dT / fluid_density;
            local_b.noalias() += eff_thermal_expansion * Tdot_int_pt * w * N;
        }
    }
}

template <typename ShapeFunction, int GlobalDim>
void StaggeredHTFEM<ShapeFunction, GlobalDim>::assembleHeatTransportEquation(
    double const t,
    double const dt,
    Eigen::VectorXd const& local_x,
    std::vector<double>& local_M_data,
    std::vector<double>& local_K_data)
{
    auto const local_p =
        local_x.template segment<pressure_size>(pressure_index);
    auto const local_T =
        local_x.template segment<temperature_size>(temperature_index);

    auto local_M = MathLib::createZeroedMatrix<LocalMatrixType>(
        local_M_data, temperature_size, temperature_size);
    auto local_K = MathLib::createZeroedMatrix<LocalMatrixType>(
        local_K_data, temperature_size, temperature_size);

    ParameterLib::SpatialPosition pos;
    pos.setElementID(this->_element.getID());

    auto const& process_data = this->_process_data;
    auto const& medium =
        *process_data.media_map.getMedium(this->_element.getID());
    auto const& liquid_phase = medium.phase("AqueousLiquid");

    auto const& b =
        process_data
            .projected_specific_body_force_vectors[this->_element.getID()];

    MaterialPropertyLib::VariableArray vars;

    unsigned const n_integration_points =
        this->_integration_method.getNumberOfPoints();

    std::vector<GlobalDimVectorType> ip_flux_vector;
    double average_velocity_norm = 0.0;
    ip_flux_vector.reserve(n_integration_points);

    for (unsigned ip(0); ip < n_integration_points; ip++)
    {
        pos.setIntegrationPoint(ip);

        auto const& ip_data = this->_ip_data[ip];
        auto const& N = ip_data.N;
        auto const& dNdx = ip_data.dNdx;
        auto const& w = ip_data.integration_weight;

        double p_at_xi = 0.;
        NumLib::shapeFunctionInterpolate(local_p, N, p_at_xi);
        double T_at_xi = 0.;
        NumLib::shapeFunctionInterpolate(local_T, N, T_at_xi);

        vars.temperature = T_at_xi;
        vars.phase_pressure = p_at_xi;

        vars.liquid_saturation = 1.0;

        auto const porosity =
            medium.property(MaterialPropertyLib::PropertyType::porosity)
                .template value<double>(vars, pos, t, dt);
        vars.porosity = porosity;

        // Use the fluid density model to compute the density
        auto const fluid_density =
            liquid_phase.property(MaterialPropertyLib::PropertyType::density)
                .template value<double>(vars, pos, t, dt);
        vars.density = fluid_density;
        auto const specific_heat_capacity_fluid =
            liquid_phase.property(MaterialPropertyLib::specific_heat_capacity)
                .template value<double>(vars, pos, t, dt);

        // Assemble mass matrix
        local_M.noalias() += w *
                             this->getHeatEnergyCoefficient(
                                 vars, porosity, fluid_density,
                                 specific_heat_capacity_fluid, pos, t, dt) *
                             N.transpose() * N;

        // Assemble Laplace matrix
        auto const viscosity =
            liquid_phase.property(MaterialPropertyLib::PropertyType::viscosity)
                .template value<double>(vars, pos, t, dt);

        auto const intrinsic_permeability =
            MaterialPropertyLib::formEigenTensor<GlobalDim>(
                medium.property(MaterialPropertyLib::PropertyType::permeability)
                    .value(vars, pos, t, dt));

        GlobalDimMatrixType const K_over_mu =
            intrinsic_permeability / viscosity;
        GlobalDimVectorType const velocity =
            process_data.has_gravity
                ? GlobalDimVectorType(-K_over_mu *
                                      (dNdx * local_p - fluid_density * b))
                : GlobalDimVectorType(-K_over_mu * dNdx * local_p);

        GlobalDimMatrixType const thermal_conductivity_dispersivity =
            this->getThermalConductivityDispersivity(
                vars, fluid_density, specific_heat_capacity_fluid, velocity,
                pos, t, dt);

        local_K.noalias() +=
            w * dNdx.transpose() * thermal_conductivity_dispersivity * dNdx;

        ip_flux_vector.emplace_back(velocity * fluid_density *
                                    specific_heat_capacity_fluid);
        average_velocity_norm += velocity.norm();
    }

    NumLib::assembleAdvectionMatrix(
        process_data.stabilizer, this->_ip_data, ip_flux_vector,
        average_velocity_norm / static_cast<double>(n_integration_points),
        local_K);
}

template <typename ShapeFunction, int GlobalDim>
std::vector<double> const&
StaggeredHTFEM<ShapeFunction, GlobalDim>::getIntPtDarcyVelocity(
    const double t,
    std::vector<GlobalVector*> const& x,
    std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
    std::vector<double>& cache) const
{
    assert(x.size() == dof_table.size());
    auto const n_processes = dof_table.size();

    std::vector<std::vector<GlobalIndexType>> indices_of_all_coupled_processes;
    indices_of_all_coupled_processes.reserve(n_processes);
    for (std::size_t process_id = 0; process_id < n_processes; ++process_id)
    {
        auto const indices =
            NumLib::getIndices(this->_element.getID(), *dof_table[process_id]);
        assert(!indices.empty());
        indices_of_all_coupled_processes.push_back(indices);
    }
    auto const local_xs =
        getCoupledLocalSolutions(x, indices_of_all_coupled_processes);

    return this->getIntPtDarcyVelocityLocal(t, local_xs, cache);
}
}  // namespace HT
}  // namespace ProcessLib
