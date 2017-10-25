/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file
 *
 * Created on August 19, 2016, 2:28 PM
 */

#pragma once

#include "LiquidFlowLocalAssembler.h"

#include "NumLib/Function/Interpolation.h"

#include "ProcessLib/StaggeredCouplingTerm.h"

#include "ProcessLib/HeatConduction/HeatConductionProcess.h"

namespace ProcessLib
{
namespace LiquidFlow
{
template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
void LiquidFlowLocalAssembler<ShapeFunction, IntegrationMethod, GlobalDim>::
    assemble(double const t, std::vector<double> const& local_x,
             std::vector<double>& local_M_data,
             std::vector<double>& local_K_data,
             std::vector<double>& local_b_data)
{
    SpatialPosition pos;
    pos.setElementID(_element.getID());
    const int material_id = _material_properties.getMaterialID(pos);

    const Eigen::MatrixXd& permeability = _material_properties.getPermeability(
        material_id, t, pos, _element.getDimension());
    // Note: For Inclined 1D in 2D/3D or 2D element in 3D, the first item in
    //  the assert must be changed to permeability.rows() ==
    //  _element->getDimension()
    assert(permeability.rows() == GlobalDim || permeability.rows() == 1);

    if (permeability.size() == 1)  // isotropic or 1D problem.
        assembleMatrixAndVector<IsotropicCalculator>(
            material_id, t, local_x, local_M_data, local_K_data, local_b_data,
            pos, permeability);
    else
        assembleMatrixAndVector<AnisotropicCalculator>(
            material_id, t, local_x, local_M_data, local_K_data, local_b_data,
            pos, permeability);
}

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
void LiquidFlowLocalAssembler<ShapeFunction, IntegrationMethod, GlobalDim>::
    assembleWithCoupledTerm(double const t, std::vector<double> const& local_x,
                            std::vector<double>& local_M_data,
                            std::vector<double>& local_K_data,
                            std::vector<double>& local_b_data,
                            LocalCouplingTerm const& coupled_term)
{
    SpatialPosition pos;
    pos.setElementID(_element.getID());
    const int material_id = _material_properties.getMaterialID(pos);

    const Eigen::MatrixXd& permeability = _material_properties.getPermeability(
        material_id, t, pos, _element.getDimension());
    // Note: For Inclined 1D in 2D/3D or 2D element in 3D, the first item in
    //  the assert must be changed to permeability.rows() ==
    //  _element->getDimension()
    assert(permeability.rows() == GlobalDim || permeability.rows() == 1);

    const double dt = coupled_term.dt;
    for (auto const& coupled_process_pair : coupled_term.coupled_processes)
    {
        if (coupled_process_pair.first ==
            std::type_index(
                typeid(ProcessLib::HeatConduction::HeatConductionProcess)))
        {
            const auto local_T0 =
                coupled_term.local_coupled_xs0.at(coupled_process_pair.first);
            const auto local_T1 =
                coupled_term.local_coupled_xs.at(coupled_process_pair.first);

            if (permeability.size() == 1)  // isotropic or 1D problem.
                assembleWithCoupledWithHeatTransport<IsotropicCalculator>(
                    material_id, t, dt, local_x, local_T0, local_T1,
                    local_M_data, local_K_data, local_b_data, pos,
                    permeability);
            else
                assembleWithCoupledWithHeatTransport<AnisotropicCalculator>(
                    material_id, t, dt, local_x, local_T0, local_T1,
                    local_M_data, local_K_data, local_b_data, pos,
                    permeability);
        }
        else
        {
            OGS_FATAL(
                "This coupled process is not presented for "
                "LiquidFlow process");
        }
    }
}

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
template <typename LaplacianGravityVelocityCalculator>
void LiquidFlowLocalAssembler<ShapeFunction, IntegrationMethod, GlobalDim>::
    assembleMatrixAndVector(const int material_id, double const t,
                            std::vector<double> const& local_x,
                            std::vector<double>& local_M_data,
                            std::vector<double>& local_K_data,
                            std::vector<double>& local_b_data,
                            SpatialPosition const& pos,
                            Eigen::MatrixXd const& permeability)
{
    auto const local_matrix_size = local_x.size();
    assert(local_matrix_size == ShapeFunction::NPOINTS);

    auto local_M = MathLib::createZeroedMatrix<NodalMatrixType>(
        local_M_data, local_matrix_size, local_matrix_size);
    auto local_K = MathLib::createZeroedMatrix<NodalMatrixType>(
        local_K_data, local_matrix_size, local_matrix_size);
    auto local_b = MathLib::createZeroedVector<NodalVectorType>(
        local_b_data, local_matrix_size);

    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    // TODO: The following two variables should be calculated inside the
    //       the integration loop for non-constant porosity and storage models.
    double porosity_variable = 0.;
    double storage_variable = 0.;
    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        auto const& sm = _shape_matrices[ip];
        auto const& wp = _integration_method.getWeightedPoint(ip);

        double p = 0.;
        NumLib::shapeFunctionInterpolate(local_x, sm.N, p);

        const double integration_factor =
            sm.integralMeasure * sm.detJ * wp.getWeight();

        // Assemble mass matrix, M
        local_M.noalias() += _material_properties.getMassCoefficient(
                                 material_id, t, pos, porosity_variable,
                                 storage_variable, p, _reference_temperature) *
                             sm.N.transpose() * sm.N * integration_factor;

        // Compute density:
        const double rho_g =
            _material_properties.getLiquidDensity(p, _reference_temperature) *
            _gravitational_acceleration;
        // Compute viscosity:
        const double mu =
            _material_properties.getViscosity(p, _reference_temperature);

        // Assemble Laplacian, K, and RHS by the gravitational term
        LaplacianGravityVelocityCalculator::calculateLaplacianAndGravityTerm(
            local_K, local_b, sm, permeability, integration_factor, mu, rho_g,
            _gravitational_axis_id);
    }
}

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
template <typename LaplacianGravityVelocityCalculator>
void LiquidFlowLocalAssembler<ShapeFunction, IntegrationMethod, GlobalDim>::
    assembleWithCoupledWithHeatTransport(
        const int material_id, double const t, double const dt,
        std::vector<double> const& local_x, std::vector<double> const& local_T0,
        std::vector<double> const& local_T1, std::vector<double>& local_M_data,
        std::vector<double>& local_K_data, std::vector<double>& local_b_data,
        SpatialPosition const& pos, Eigen::MatrixXd const& permeability)
{
    auto const local_matrix_size = local_x.size();
    assert(local_matrix_size == ShapeFunction::NPOINTS);

    auto local_M = MathLib::createZeroedMatrix<NodalMatrixType>(
        local_M_data, local_matrix_size, local_matrix_size);
    auto local_K = MathLib::createZeroedMatrix<NodalMatrixType>(
        local_K_data, local_matrix_size, local_matrix_size);
    auto local_b = MathLib::createZeroedVector<NodalVectorType>(
        local_b_data, local_matrix_size);

    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    // TODO: The following two variables should be calculated inside the
    //       the integration loop for non-constant porosity and storage models.
    double porosity_variable = 0.;
    double storage_variable = 0.;
    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        auto const& sm = _shape_matrices[ip];
        auto const& wp = _integration_method.getWeightedPoint(ip);

        double p = 0.;
        NumLib::shapeFunctionInterpolate(local_x, sm.N, p);
        double T0 = 0.;
        NumLib::shapeFunctionInterpolate(local_T0, sm.N, T0);
        double T = 0.;
        NumLib::shapeFunctionInterpolate(local_T1, sm.N, T);

        const double integration_factor =
            sm.integralMeasure * sm.detJ * wp.getWeight();

        // Assemble mass matrix, M
        local_M.noalias() += _material_properties.getMassCoefficient(
                                 material_id, t, pos, porosity_variable,
                                 storage_variable, p, T) *
                             sm.N.transpose() * sm.N * integration_factor;

        // Compute density:
        const double rho = _material_properties.getLiquidDensity(p, T);
        const double rho_g = rho * _gravitational_acceleration;

        // Compute viscosity:
        const double mu = _material_properties.getViscosity(p, T);

        // Assemble Laplacian, K, and RHS by the gravitational term
        LaplacianGravityVelocityCalculator::calculateLaplacianAndGravityTerm(
            local_K, local_b, sm, permeability, integration_factor, mu, rho_g,
            _gravitational_axis_id);

        // Add the thermal expansion term
        auto const solid_thermal_expansion =
            _material_properties.getSolidThermalExpansion(t, pos);
        auto const biot_constant = _material_properties.getBiotConstant(t, pos);
        auto const porosity = _material_properties.getPorosity(
            material_id, t, pos, porosity_variable, T);
        const double eff_thermal_expansion =
            3.0 * (biot_constant - porosity) * solid_thermal_expansion -
            porosity * _material_properties.getdLiquidDensity_dT(p, T) / rho;
        local_b.noalias() +=
            eff_thermal_expansion * (T - T0) * integration_factor * sm.N / dt;
    }
}

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
std::vector<double> const&
LiquidFlowLocalAssembler<ShapeFunction, IntegrationMethod, GlobalDim>::
    getIntPtDarcyVelocity(const double t,
                          GlobalVector const& current_solution,
                          NumLib::LocalToGlobalIndexMap const& dof_table,
                          std::vector<double>& veloctiy_cache) const
{
    auto const indices = NumLib::getIndices(_element.getID(), dof_table);
    auto const local_x = current_solution.get(indices);
    auto const num_intpts = _integration_method.getNumberOfPoints();
    veloctiy_cache.clear();
    auto veloctiy_cache_vectors = MathLib::createZeroedMatrix<
        Eigen::Matrix<double, GlobalDim, Eigen::Dynamic, Eigen::RowMajor>>(
        veloctiy_cache, GlobalDim, num_intpts);

    SpatialPosition pos;
    pos.setElementID(_element.getID());
    const int material_id = _material_properties.getMaterialID(pos);
    const Eigen::MatrixXd& permeability = _material_properties.getPermeability(
        material_id, t, pos, _element.getDimension());

    // Note: For Inclined 1D in 2D/3D or 2D element in 3D, the first item in
    //  the assert must be changed to perm.rows() == _element->getDimension()
    assert(permeability.rows() == GlobalDim || permeability.rows() == 1);

    if (!_coupling_term)
    {
        computeDarcyVelocity(permeability, local_x, veloctiy_cache_vectors);
    }
    else
    {
        auto const local_coupled_xs =
            getCurrentLocalSolutionsOfCoupledProcesses(
                _coupling_term->coupled_xs, indices);
        if (local_coupled_xs.empty())
            computeDarcyVelocity(permeability, local_x, veloctiy_cache_vectors);
        else
            computeDarcyVelocityWithCoupling(permeability, local_x,
                                             local_coupled_xs,
                                             veloctiy_cache_vectors);
    }

    return veloctiy_cache;
}

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
void LiquidFlowLocalAssembler<ShapeFunction, IntegrationMethod, GlobalDim>::
    computeDarcyVelocity(
        Eigen::MatrixXd const& permeability,
        std::vector<double> const& local_x,
        MatrixOfVelocityAtIntegrationPoints& darcy_velocity_at_ips) const
{
    if (permeability.size() == 1)  // isotropic or 1D problem.
        computeDarcyVelocityLocal<IsotropicCalculator>(local_x, permeability,
                                                       darcy_velocity_at_ips);
    else
        computeDarcyVelocityLocal<AnisotropicCalculator>(local_x, permeability,
                                                         darcy_velocity_at_ips);
}

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
template <typename LaplacianGravityVelocityCalculator>
void LiquidFlowLocalAssembler<ShapeFunction, IntegrationMethod, GlobalDim>::
    computeDarcyVelocityLocal(
        std::vector<double> const& local_x,
        Eigen::MatrixXd const& permeability,
        MatrixOfVelocityAtIntegrationPoints& darcy_velocity_at_ips) const
{
    auto const local_matrix_size = local_x.size();
    assert(local_matrix_size == ShapeFunction::NPOINTS);

    const auto local_p_vec =
        MathLib::toVector<NodalVectorType>(local_x, local_matrix_size);

    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        auto const& sm = _shape_matrices[ip];
        double p = 0.;
        NumLib::shapeFunctionInterpolate(local_x, sm.N, p);

        const double rho_g =
            _material_properties.getLiquidDensity(p, _reference_temperature) *
            _gravitational_acceleration;
        // Compute viscosity:
        const double mu =
            _material_properties.getViscosity(p, _reference_temperature);

        LaplacianGravityVelocityCalculator::calculateVelocity(
            ip, local_p_vec, sm, permeability, mu, rho_g,
            _gravitational_axis_id, darcy_velocity_at_ips);
    }
}

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
void LiquidFlowLocalAssembler<ShapeFunction, IntegrationMethod, GlobalDim>::
    computeDarcyVelocityWithCoupling(
        Eigen::MatrixXd const& permeability, std::vector<double> const& local_x,
        std::unordered_map<std::type_index, const std::vector<double>> const&
            coupled_local_solutions,
        MatrixOfVelocityAtIntegrationPoints& darcy_velocity_at_ips) const
{
    const auto local_T = coupled_local_solutions.at(std::type_index(
        typeid(ProcessLib::HeatConduction::HeatConductionProcess)));

    if (permeability.size() == 1)  // isotropic or 1D problem.
        computeDarcyVelocityCoupledWithHeatTransportLocal<IsotropicCalculator>(
            local_x, local_T, permeability, darcy_velocity_at_ips);
    else
        computeDarcyVelocityCoupledWithHeatTransportLocal<
            AnisotropicCalculator>(local_x, local_T, permeability,
                                   darcy_velocity_at_ips);
}

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
template <typename LaplacianGravityVelocityCalculator>
void LiquidFlowLocalAssembler<ShapeFunction, IntegrationMethod, GlobalDim>::
    computeDarcyVelocityCoupledWithHeatTransportLocal(
        std::vector<double> const& local_x,
        std::vector<double> const& local_T,
        Eigen::MatrixXd const& permeability,
        MatrixOfVelocityAtIntegrationPoints& darcy_velocity_at_ips) const
{
    auto const local_matrix_size = local_x.size();
    assert(local_matrix_size == ShapeFunction::NPOINTS);

    const auto local_p_vec =
        MathLib::toVector<NodalVectorType>(local_x, local_matrix_size);

    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        auto const& sm = _shape_matrices[ip];
        double p = 0.;
        NumLib::shapeFunctionInterpolate(local_x, sm.N, p);
        double T = 0.;
        NumLib::shapeFunctionInterpolate(local_T, sm.N, T);

        const double rho_g = _material_properties.getLiquidDensity(p, T) *
                             _gravitational_acceleration;
        // Compute viscosity:
        const double mu = _material_properties.getViscosity(p, T);

        LaplacianGravityVelocityCalculator::calculateVelocity(
            ip, local_p_vec, sm, permeability, mu, rho_g,
            _gravitational_axis_id, darcy_velocity_at_ips);
    }
}

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
void LiquidFlowLocalAssembler<ShapeFunction, IntegrationMethod, GlobalDim>::
    IsotropicCalculator::calculateLaplacianAndGravityTerm(
        Eigen::Map<NodalMatrixType>& local_K,
        Eigen::Map<NodalVectorType>& local_b, ShapeMatrices const& sm,
        Eigen::MatrixXd const& permeability, double const integration_factor,
        double const mu, double const rho_g, int const gravitational_axis_id)
{
    const double K = permeability(0, 0) / mu;
    const double fac = K * integration_factor;
    local_K.noalias() += fac * sm.dNdx.transpose() * sm.dNdx;

    if (gravitational_axis_id >= 0)
    {
        // Note: Since permeability, K, is a scalar number in this function,
        // (gradN)*K*V becomes K*(grad N)*V. With the gravity vector of V,
        // the simplification of (grad N)*V is the gravitational_axis_id th
        // column of the transpose of (grad N) multiplied with rho_g.
        local_b.noalias() -=
            fac * sm.dNdx.transpose().col(gravitational_axis_id) * rho_g;
    }
}

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
void LiquidFlowLocalAssembler<ShapeFunction, IntegrationMethod, GlobalDim>::
    IsotropicCalculator::calculateVelocity(
        unsigned const ip, Eigen::Map<const NodalVectorType> const& local_p,
        ShapeMatrices const& sm, Eigen::MatrixXd const& permeability,
        double const mu, double const rho_g, int const gravitational_axis_id,
        MatrixOfVelocityAtIntegrationPoints& darcy_velocity_at_ips)
{
    const double K = permeability(0, 0) / mu;
    // Compute the velocity
    darcy_velocity_at_ips.col(ip).noalias() = -K * sm.dNdx * local_p;
    // gravity term
    if (gravitational_axis_id >= 0)
        darcy_velocity_at_ips.col(ip)[gravitational_axis_id] -= K * rho_g;
}

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
void LiquidFlowLocalAssembler<ShapeFunction, IntegrationMethod, GlobalDim>::
    AnisotropicCalculator::calculateLaplacianAndGravityTerm(
        Eigen::Map<NodalMatrixType>& local_K,
        Eigen::Map<NodalVectorType>& local_b, ShapeMatrices const& sm,
        Eigen::MatrixXd const& permeability, double const integration_factor,
        double const mu, double const rho_g, int const gravitational_axis_id)
{
    const double fac = integration_factor / mu;
    local_K.noalias() += fac * sm.dNdx.transpose() * permeability * sm.dNdx;
    if (gravitational_axis_id >= 0)
    {
        // Note: permeability * gravity_vector = rho_g *
        // permeability.col(gravitational_axis_id)
        //       i.e the result is the gravitational_axis_id th column of
        //       permeability multiplied with rho_g.
        local_b.noalias() -= fac * rho_g * sm.dNdx.transpose() *
                             permeability.col(gravitational_axis_id);
    }
}

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
void LiquidFlowLocalAssembler<ShapeFunction, IntegrationMethod, GlobalDim>::
    AnisotropicCalculator::calculateVelocity(
        unsigned const ip, Eigen::Map<const NodalVectorType> const& local_p,
        ShapeMatrices const& sm, Eigen::MatrixXd const& permeability,
        double const mu, double const rho_g, int const gravitational_axis_id,
        MatrixOfVelocityAtIntegrationPoints& darcy_velocity_at_ips)
{
    // Compute the velocity
    darcy_velocity_at_ips.col(ip).noalias() =
        -permeability * sm.dNdx * local_p / mu;
    if (gravitational_axis_id >= 0)
    {
        darcy_velocity_at_ips.col(ip).noalias() -=
            rho_g * permeability.col(gravitational_axis_id) / mu;
    }
}

}  // end of namespace
}  // end of namespace
