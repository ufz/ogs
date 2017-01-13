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

#include "MaterialLib/PhysicalConstant.h"
#include "NumLib/Function/Interpolation.h"

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

    const Eigen::MatrixXd& perm = _material_properties.getPermeability(
        material_id, t, pos, _element.getDimension());
    // Note: For Inclined 1D in 2D/3D or 2D element in 3D, the first item in
    //  the assert must be changed to perm.rows() == _element->getDimension()
    assert(perm.rows() == GlobalDim || perm.rows() == 1);

    if (perm.size() == 1)  // isotropic or 1D problem.
        local_assemble<IsotropicCalculator>(material_id, t, local_x,
                                            local_M_data, local_K_data,
                                            local_b_data, pos, perm);
    else
        local_assemble<AnisotropicCalculator>(material_id, t, local_x,
                                              local_M_data, local_K_data,
                                              local_b_data, pos, perm);
}

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
void LiquidFlowLocalAssembler<ShapeFunction, IntegrationMethod, GlobalDim>::
    coupling_assemble(double const t, std::vector<double> const& local_x,
                      std::vector<double>& local_M_data,
                      std::vector<double>& local_K_data,
                      std::vector<double>& local_b_data,
                      LocalCouplingTerm const& coupled_term)
{
    SpatialPosition pos;
    pos.setElementID(_element.getID());
    const int material_id = _material_properties.getMaterialID(pos);

    const Eigen::MatrixXd& perm = _material_properties.getPermeability(
        material_id, t, pos, _element.getDimension());
    // Note: For Inclined 1D in 2D/3D or 2D element in 3D, the first item in
    //  the assert must be changed to perm.rows() == _element->getDimension()
    assert(perm.rows() == GlobalDim || perm.rows() == 1);

    const double dt = coupled_term.dt;
    auto it = coupled_term.coupled_processes.begin();
    while (it != coupled_term.coupled_processes.end())
    {
        switch (it->first)
        {
            case ProcessLib::ProcessType::HeatConductionProcess:
            {
                const auto local_T0 = coupled_term.local_coupled_xs0.at(
                    ProcessLib::ProcessType::HeatConductionProcess);
                const auto local_T1 = coupled_term.local_coupled_xs.at(
                    ProcessLib::ProcessType::HeatConductionProcess);

                if (perm.size() == 1)  // isotropic or 1D problem.
                    local_assembleCoupledWithHeatTransport<IsotropicCalculator>(
                        material_id, t, dt, local_x, local_T0, local_T1,
                        local_M_data, local_K_data, local_b_data, pos, perm);
                else
                    local_assembleCoupledWithHeatTransport<
                        AnisotropicCalculator>(
                        material_id, t, dt, local_x, local_T0, local_T1,
                        local_M_data, local_K_data, local_b_data, pos, perm);
            }
            break;
            default:
                OGS_FATAL(
                    "This coupled process is not presented for "
                    "LiquidFlow process");
        }
        it++;
    }
}

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
template <typename LaplacianGravityVelocityCalculator>
void LiquidFlowLocalAssembler<ShapeFunction, IntegrationMethod, GlobalDim>::
    local_assemble(const int material_id, double const t,
                   std::vector<double> const& local_x,
                   std::vector<double>& local_M_data,
                   std::vector<double>& local_K_data,
                   std::vector<double>& local_b_data,
                   SpatialPosition const& pos, Eigen::MatrixXd const& perm)
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
    const double temperature =
        MaterialLib::PhysicalConstant::CelsiusZeroInKelvin + 18.0;
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
                                 storage_variable, p, temperature) *
                             sm.N.transpose() * sm.N * integration_factor;

        // Compute density:
        const double rho_g =
            _material_properties.getLiquidDensity(p, temperature) *
            _gravitational_acceleration;
        // Compute viscosity:
        const double mu = _material_properties.getViscosity(p, temperature);

        // Assemble Laplacian, K, and RHS by the gravitational term
        LaplacianGravityVelocityCalculator::calculateLaplacianAndGravityTerm(
            local_K, local_b, sm, perm, integration_factor, mu, rho_g,
            _gravitational_axis_id);
    }
}

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
template <typename LaplacianGravityVelocityCalculator>
void LiquidFlowLocalAssembler<ShapeFunction, IntegrationMethod, GlobalDim>::
    local_assembleCoupledWithHeatTransport(
        const int material_id, double const t, double const dt,
        std::vector<double> const& local_x, std::vector<double> const& local_T0,
        std::vector<double> const& local_T1, std::vector<double>& local_M_data,
        std::vector<double>& local_K_data, std::vector<double>& local_b_data,
        SpatialPosition const& pos, Eigen::MatrixXd const& perm)
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
            local_K, local_b, sm, perm, integration_factor, mu, rho_g,
            _gravitational_axis_id);

        // Add the thermal expansion term
        auto const solid_thermal_expansion =
            _material_properties.getSolidThermalExpansion(t, pos);
        auto const biot_constant =
            _material_properties.getBiotConstant(t, pos);
        auto const porosity = _material_properties.getPorosity(
            material_id, porosity_variable, T);
        const double eff_thermal_expansion =
            3.0 * (biot_constant - porosity) * solid_thermal_expansion +
            porosity * _material_properties.getdLiquidDensity_dT(p, T) / rho;
        local_b.noalias() +=
            eff_thermal_expansion * (T - T0) * integration_factor * sm.N / dt;
    }
}

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
void LiquidFlowLocalAssembler<ShapeFunction, IntegrationMethod, GlobalDim>::
    computeSecondaryVariableConcrete(double const t,
                                     std::vector<double> const& local_x)
{
    SpatialPosition pos;
    pos.setElementID(_element.getID());
    const int material_id = _material_properties.getMaterialID(pos);
    const Eigen::MatrixXd& perm = _material_properties.getPermeability(
        material_id, t, pos, _element.getDimension());

    // Note: For Inclined 1D in 2D/3D or 2D element in 3D, the first item in
    //  the assert must be changed to perm.rows() == _element->getDimension()
    assert(perm.rows() == GlobalDim || perm.rows() == 1);

    if (perm.size() == 1)  // isotropic or 1D problem.
        computeSecondaryVariableLocal<IsotropicCalculator>(t, local_x, pos,
                                                           perm);
    else
        computeSecondaryVariableLocal<AnisotropicCalculator>(t, local_x, pos,
                                                             perm);
}

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
template <typename LaplacianGravityVelocityCalculator>
void LiquidFlowLocalAssembler<ShapeFunction, IntegrationMethod, GlobalDim>::
    computeSecondaryVariableLocal(double const /*t*/,
                                  std::vector<double> const& local_x,
                                  SpatialPosition const& /*pos*/,
                                  Eigen::MatrixXd const& perm)
{
    auto const local_matrix_size = local_x.size();
    assert(local_matrix_size == ShapeFunction::NPOINTS);

    const auto local_p_vec =
        MathLib::toVector<NodalVectorType>(local_x, local_matrix_size);

    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    const double temperature =
        MaterialLib::PhysicalConstant::CelsiusZeroInKelvin + 18.0;
    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        auto const& sm = _shape_matrices[ip];
        double p = 0.;
        NumLib::shapeFunctionInterpolate(local_x, sm.N, p);

        const double rho_g =
            _material_properties.getLiquidDensity(p, temperature) *
            _gravitational_acceleration;
        // Compute viscosity:
        const double mu = _material_properties.getViscosity(p, temperature);

        LaplacianGravityVelocityCalculator::calculateVelocity(
            _darcy_velocities, local_p_vec, sm, perm, ip, mu, rho_g,
            _gravitational_axis_id);
    }
}

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
template <typename LaplacianGravityVelocityCalculator>
void LiquidFlowLocalAssembler<ShapeFunction, IntegrationMethod, GlobalDim>::
    computeSecondaryVariableCoupledWithHeatTransportLocal(
        double const /*t*/,
        std::vector<double> const& local_x,
        std::vector<double> const& local_T,
        SpatialPosition const& /*pos*/,
        Eigen::MatrixXd const& perm)
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
            _darcy_velocities, local_p_vec, sm, perm, ip, mu, rho_g,
            _gravitational_axis_id);
    }
}

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
void LiquidFlowLocalAssembler<ShapeFunction, IntegrationMethod, GlobalDim>::
    IsotropicCalculator::calculateLaplacianAndGravityTerm(
        Eigen::Map<NodalMatrixType>& local_K,
        Eigen::Map<NodalVectorType>& local_b, ShapeMatrices const& sm,
        Eigen::MatrixXd const& perm, double const integration_factor,
        double const mu, double const rho_g, int const gravitational_axis_id)
{
    const double K = perm(0, 0) / mu;
    const double fac = K * integration_factor;
    local_K.noalias() += fac * sm.dNdx.transpose() * sm.dNdx;

    if (gravitational_axis_id >= 0)
    {
        // Note: Since perm, K, is a scalar number in this function,
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
        std::vector<std::vector<double>>& darcy_velocities,
        Eigen::Map<const NodalVectorType> const& local_p,
        ShapeMatrices const& sm, Eigen::MatrixXd const& perm, unsigned const ip,
        double const mu, double const rho_g, int const gravitational_axis_id)
{
    const double K = perm(0, 0) / mu;
    // Compute the velocity
    GlobalDimVectorType darcy_velocity = -K * sm.dNdx * local_p;
    // gravity term
    if (gravitational_axis_id >= 0)
        darcy_velocity[gravitational_axis_id] -= K * rho_g;
    for (unsigned d = 0; d < GlobalDim; ++d)
    {
        darcy_velocities[d][ip] = darcy_velocity[d];
    }
}

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
void LiquidFlowLocalAssembler<ShapeFunction, IntegrationMethod, GlobalDim>::
    AnisotropicCalculator::calculateLaplacianAndGravityTerm(
        Eigen::Map<NodalMatrixType>& local_K,
        Eigen::Map<NodalVectorType>& local_b, ShapeMatrices const& sm,
        Eigen::MatrixXd const& perm, double const integration_factor,
        double const mu, double const rho_g, int const gravitational_axis_id)
{
    const double fac = integration_factor / mu;
    local_K.noalias() += fac * sm.dNdx.transpose() * perm * sm.dNdx;
    if (gravitational_axis_id >= 0)
    {
        // Note: perm * gravity_vector = rho_g * perm.col(gravitational_axis_id)
        //       i.e the result is the gravitational_axis_id th column of
        //       perm multiplied with rho_g.
        local_b.noalias() -=
            fac * rho_g * sm.dNdx.transpose() * perm.col(gravitational_axis_id);
    }
}

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
void LiquidFlowLocalAssembler<ShapeFunction, IntegrationMethod, GlobalDim>::
    AnisotropicCalculator::calculateVelocity(
        std::vector<std::vector<double>>& darcy_velocities,
        Eigen::Map<const NodalVectorType> const& local_p,
        ShapeMatrices const& sm, Eigen::MatrixXd const& perm, unsigned const ip,
        double const mu, double const rho_g, int const gravitational_axis_id)
{
    // Compute the velocity
    GlobalDimVectorType darcy_velocity = -perm * sm.dNdx * local_p / mu;
    if (gravitational_axis_id >= 0)
    {
        darcy_velocity.noalias() -=
            rho_g * perm.col(gravitational_axis_id) / mu;
    }
    for (unsigned d = 0; d < GlobalDim; ++d)
    {
        darcy_velocities[d][ip] = darcy_velocity[d];
    }
}

}  // end of namespace
}  // end of namespace
