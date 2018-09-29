/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  \file   StaggeredHTFEM-impl.h
 *  Created on October 13, 2017, 3:52 PM
 */

#pragma once

#include "StaggeredHTFEM.h"

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

    SpatialPosition pos;
    pos.setElementID(this->_element.getID());

    auto const& material_properties = this->_material_properties;

    auto const& b = material_properties.specific_body_force;

    MaterialLib::Fluid::FluidProperty::ArrayType vars;

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

        auto const porosity =
            material_properties.porous_media_properties.getPorosity(t, pos)
                .getValue(t, pos, 0.0, T1_at_xi);

        vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::T)] =
            T1_at_xi;
        vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::p)] =
            p_at_xi;

        // Use the fluid density model to compute the density
        auto const fluid_density =
            material_properties.fluid_properties->getValue(
                MaterialLib::Fluid::FluidPropertyType::Density, vars);
        const double dfluid_density_dp =
            material_properties.fluid_properties->getdValue(
                MaterialLib::Fluid::FluidPropertyType::Density, vars,
                MaterialLib::Fluid::PropertyVariableType::p);

        // Use the viscosity model to compute the viscosity
        auto const viscosity = material_properties.fluid_properties->getValue(
            MaterialLib::Fluid::FluidPropertyType::Viscosity, vars);

        // \todo the argument to getValue() has to be changed for non
        // constant storage model
        auto const specific_storage =
            material_properties.porous_media_properties
                .getSpecificStorage(t, pos)
                .getValue(0.0);

        auto const intrinsic_permeability =
            material_properties.porous_media_properties
                .getIntrinsicPermeability(t, pos)
                .getValue(t, pos, 0.0, T1_at_xi);
        GlobalDimMatrixType K_over_mu = intrinsic_permeability / viscosity;

        // matrix assembly
        local_M.noalias() +=
            w *
            (porosity * dfluid_density_dp / fluid_density + specific_storage) *
            N.transpose() * N;

        local_K.noalias() += w * dNdx.transpose() * K_over_mu * dNdx;

        if (material_properties.has_gravity)
        {
            local_b.noalias() +=
                w * fluid_density * dNdx.transpose() * K_over_mu * b;
        }

        if (!material_properties.has_fluid_thermal_expansion)
            return;

        // Add the thermal expansion term
        {
            auto const solid_thermal_expansion =
                material_properties.solid_thermal_expansion(t, pos)[0];
            const double dfluid_density_dT =
                material_properties.fluid_properties->getdValue(
                    MaterialLib::Fluid::FluidPropertyType::Density, vars,
                    MaterialLib::Fluid::PropertyVariableType::T);
            double T0_at_xi = 0.;
            NumLib::shapeFunctionInterpolate(local_T0, N, T0_at_xi);
            auto const biot_constant =
                material_properties.biot_constant(t, pos)[0];
            const double eff_thermal_expansion =
                3.0 * (biot_constant - porosity) * solid_thermal_expansion -
                porosity * dfluid_density_dT / fluid_density;
            local_b.noalias() +=
                eff_thermal_expansion * (T1_at_xi - T0_at_xi) * w * N / dt;
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

    SpatialPosition pos;
    pos.setElementID(this->_element.getID());

    auto const& material_properties = this->_material_properties;

    auto const& b = material_properties.specific_body_force;

    GlobalDimMatrixType const& I(
        GlobalDimMatrixType::Identity(GlobalDim, GlobalDim));

    MaterialLib::Fluid::FluidProperty::ArrayType vars;

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

        auto const porosity =
            material_properties.porous_media_properties.getPorosity(t, pos)
                .getValue(t, pos, 0.0, T1_at_xi);

        vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::T)] =
            T1_at_xi;
        vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::p)] =
            p_at_xi;
        auto const fluid_density =
            material_properties.fluid_properties->getValue(
                MaterialLib::Fluid::FluidPropertyType::Density, vars);
        auto const specific_heat_capacity_fluid =
            material_properties.fluid_properties->getValue(
                MaterialLib::Fluid::FluidPropertyType::HeatCapacity, vars);

        // Assemble mass matrix
        local_M.noalias() +=
            w *
            this->getHeatEnergyCoefficient(t, pos, porosity, fluid_density,
                                           specific_heat_capacity_fluid) *
            N.transpose() * N;

        // Assemble Laplace matrix
        auto const viscosity = material_properties.fluid_properties->getValue(
            MaterialLib::Fluid::FluidPropertyType::Viscosity, vars);

        auto const intrinsic_permeability =
            material_properties.porous_media_properties
                .getIntrinsicPermeability(t, pos)
                .getValue(t, pos, 0.0, T1_at_xi);

        GlobalDimMatrixType K_over_mu = intrinsic_permeability / viscosity;
        GlobalDimVectorType const velocity =
            material_properties.has_gravity
                ? GlobalDimVectorType(-K_over_mu * (dNdx * local_p_Eigen_type -
                                                    fluid_density * b))
                : GlobalDimVectorType(-K_over_mu * dNdx * local_p_Eigen_type);

        GlobalDimMatrixType const thermal_conductivity_dispersivity =
            this->getThermalConductivityDispersivity(
                t, pos, porosity, fluid_density, specific_heat_capacity_fluid,
                velocity, I);

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
