/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 * Created on August 19, 2016, 2:28 PM
 */

#pragma once

#include "LiquidFlowLocalAssembler.h"
#include "MaterialLib/MPL/MaterialSpatialDistributionMap.h"
#include "MaterialLib/MPL/Utils/FormEigenTensor.h"
#include "MaterialLib/MPL/Utils/GetFluidDensityAndViscosity.h"
#include "MaterialLib/MPL/VariableType.h"
#include "NumLib/Fem/Interpolation.h"

namespace ProcessLib
{
namespace LiquidFlow
{
template <typename ShapeFunction, int GlobalDim>
void LiquidFlowLocalAssembler<ShapeFunction, GlobalDim>::assemble(
    double const t, double const dt, std::vector<double> const& local_x,
    std::vector<double> const& /*local_x_prev*/,
    std::vector<double>& local_M_data, std::vector<double>& local_K_data,
    std::vector<double>& local_b_data)
{
    ParameterLib::SpatialPosition pos;
    pos.setElementID(_element.getID());

    auto const& medium = *_process_data.media_map.getMedium(_element.getID());
    MaterialPropertyLib::VariableArray vars;
    vars.temperature =
        medium[MaterialPropertyLib::PropertyType::reference_temperature]
            .template value<double>(vars, pos, t, dt);
    vars.liquid_phase_pressure = std::numeric_limits<double>::quiet_NaN();
    GlobalDimMatrixType const permeability =
        MaterialPropertyLib::formEigenTensor<GlobalDim>(
            medium[MaterialPropertyLib::PropertyType::permeability].value(
                vars, pos, t, dt));
    // Note: For Inclined 1D in 2D/3D or 2D element in 3D, the first item in
    //  the assert must be changed to permeability.rows() ==
    //  _element->getDimension()
    assert(permeability.rows() == GlobalDim || permeability.rows() == 1);

    if (permeability.size() == 1)
    {  // isotropic or 1D problem.
        assembleMatrixAndVector<IsotropicCalculator>(
            t, dt, local_x, local_M_data, local_K_data, local_b_data);
    }
    else
    {
        assembleMatrixAndVector<AnisotropicCalculator>(
            t, dt, local_x, local_M_data, local_K_data, local_b_data);
    }
}

template <typename ShapeFunction, int GlobalDim>
Eigen::Vector3d LiquidFlowLocalAssembler<ShapeFunction, GlobalDim>::getFlux(
    MathLib::Point3d const& p_local_coords, double const t,
    std::vector<double> const& local_x) const
{
    // TODO (tf) Temporary value not used by current material models. Need
    // extension of getFlux interface
    double const dt = std::numeric_limits<double>::quiet_NaN();

    // Note: Axial symmetry is set to false here, because we only need dNdx
    // here, which is not affected by axial symmetry.
    auto const shape_matrices =
        NumLib::computeShapeMatrices<ShapeFunction, ShapeMatricesType,
                                     GlobalDim>(_element,
                                                false /*is_axially_symmetric*/,
                                                std::array{p_local_coords})[0];

    // create pos object to access the correct media property
    ParameterLib::SpatialPosition pos;
    pos.setElementID(_element.getID());

    auto const& medium = *_process_data.media_map.getMedium(_element.getID());
    auto const& fluid_phase = fluidPhase(medium);

    MaterialPropertyLib::VariableArray vars;

    double pressure = 0.0;
    NumLib::shapeFunctionInterpolate(local_x, shape_matrices.N, pressure);
    vars.liquid_phase_pressure = pressure;

    GlobalDimMatrixType const intrinsic_permeability =
        MaterialPropertyLib::formEigenTensor<GlobalDim>(
            medium[MaterialPropertyLib::PropertyType::permeability].value(
                vars, pos, t, dt));
    auto const viscosity =
        fluid_phase[MaterialPropertyLib::PropertyType::viscosity]
            .template value<double>(vars, pos, t, dt);

    Eigen::Vector3d flux(0.0, 0.0, 0.0);
    flux.head<GlobalDim>() =
        -intrinsic_permeability / viscosity * shape_matrices.dNdx *
        Eigen::Map<const NodalVectorType>(local_x.data(), local_x.size());

    return flux;
}

template <typename ShapeFunction, int GlobalDim>
template <typename LaplacianGravityVelocityCalculator>
void LiquidFlowLocalAssembler<ShapeFunction, GlobalDim>::
    assembleMatrixAndVector(double const t, double const dt,
                            std::vector<double> const& local_x,
                            std::vector<double>& local_M_data,
                            std::vector<double>& local_K_data,
                            std::vector<double>& local_b_data)
{
    auto const local_matrix_size = local_x.size();
    assert(local_matrix_size == ShapeFunction::NPOINTS);

    auto local_M = MathLib::createZeroedMatrix<NodalMatrixType>(
        local_M_data, local_matrix_size, local_matrix_size);
    auto local_K = MathLib::createZeroedMatrix<NodalMatrixType>(
        local_K_data, local_matrix_size, local_matrix_size);
    auto local_b = MathLib::createZeroedVector<NodalVectorType>(
        local_b_data, local_matrix_size);
    const auto local_p_vec =
        MathLib::toVector<NodalVectorType>(local_x, local_matrix_size);

    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    ParameterLib::SpatialPosition pos;
    pos.setElementID(_element.getID());

    auto const& medium = *_process_data.media_map.getMedium(_element.getID());
    auto const& fluid_phase = fluidPhase(medium);

    MaterialPropertyLib::VariableArray vars;
    auto& phase_pressure = medium.hasPhase("Gas") ? vars.gas_phase_pressure
                                                  : vars.liquid_phase_pressure;

    vars.temperature =
        medium[MaterialPropertyLib::PropertyType::reference_temperature]
            .template value<double>(vars, pos, t, dt);

    GlobalDimVectorType const projected_body_force_vector =
        _process_data.element_rotation_matrices[_element.getID()] *
        _process_data.element_rotation_matrices[_element.getID()].transpose() *
        _process_data.specific_body_force;

    auto const& Ns = _process_data.shape_matrix_cache
                         .NsHigherOrder<typename ShapeFunction::MeshElement>();

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        auto const& ip_data = _ip_data[ip];
        auto const& N = Ns[ip];

        phase_pressure = N.dot(local_p_vec);

        auto const [fluid_density, viscosity] =
            MaterialPropertyLib::getFluidDensityAndViscosity(t, dt, pos,
                                                             fluid_phase, vars);

        auto const porosity =
            medium[MaterialPropertyLib::PropertyType::porosity]
                .template value<double>(vars, pos, t, dt);
        auto const specific_storage =
            medium[MaterialPropertyLib::PropertyType::storage]
                .template value<double>(vars, pos, t, dt);

        auto const get_drho_dp = [&]()
        {
            return fluid_phase[MaterialPropertyLib::PropertyType::density]
                .template dValue<double>(vars, _process_data.phase_variable,
                                         pos, t, dt);
        };
        auto const storage =
            _process_data.equation_balance_type == EquationBalanceType::volume
                ? specific_storage
                : specific_storage + porosity * get_drho_dp() / fluid_density;

        double const scaling_factor =
            _process_data.equation_balance_type == EquationBalanceType::volume
                ? 1.0
                : fluid_density;
        // Assemble mass matrix, M
        local_M.noalias() += scaling_factor * storage * N.transpose() * N *
                             ip_data.integration_weight;

        GlobalDimMatrixType const permeability =
            MaterialPropertyLib::formEigenTensor<GlobalDim>(
                medium[MaterialPropertyLib::PropertyType::permeability].value(
                    vars, pos, t, dt));

        // Assemble Laplacian, K, and RHS by the gravitational term
        LaplacianGravityVelocityCalculator::calculateLaplacianAndGravityTerm(
            local_K, local_b, ip_data, scaling_factor * permeability, viscosity,
            fluid_density, projected_body_force_vector,
            _process_data.has_gravity);
    }
}

template <typename ShapeFunction, int GlobalDim>
template <typename VelocityCacheType>
void LiquidFlowLocalAssembler<ShapeFunction, GlobalDim>::computeDarcyVelocity(
    bool const is_scalar_permeability, const double t, const double dt,
    std::vector<double> const& local_x,
    ParameterLib::SpatialPosition const& pos,
    VelocityCacheType& darcy_velocity_at_ips) const
{
    if (is_scalar_permeability)
    {  // isotropic or 1D problem.
        computeProjectedDarcyVelocity<IsotropicCalculator>(
            t, dt, local_x, pos, darcy_velocity_at_ips);
    }
    else
    {
        computeProjectedDarcyVelocity<AnisotropicCalculator>(
            t, dt, local_x, pos, darcy_velocity_at_ips);
    }
}

template <typename ShapeFunction, int GlobalDim>
std::vector<double> const&
LiquidFlowLocalAssembler<ShapeFunction, GlobalDim>::getIntPtDarcyVelocity(
    const double t,
    std::vector<GlobalVector*> const& x,
    std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_tables,
    std::vector<double>& velocity_cache) const
{
    // TODO (tf) Temporary value not used by current material models. Need
    // extension of secondary variable interface.
    double const dt = std::numeric_limits<double>::quiet_NaN();

    constexpr int process_id = 0;
    auto const indices =
        NumLib::getIndices(_element.getID(), *dof_tables[process_id]);
    auto const local_x = x[process_id]->get(indices);
    auto const n_integration_points = _integration_method.getNumberOfPoints();
    velocity_cache.clear();

    ParameterLib::SpatialPosition pos;
    pos.setElementID(_element.getID());

    auto const& medium = *_process_data.media_map.getMedium(_element.getID());
    MaterialPropertyLib::VariableArray vars;
    vars.temperature =
        medium[MaterialPropertyLib::PropertyType::reference_temperature]
            .template value<double>(vars, pos, t, dt);
    vars.liquid_phase_pressure = std::numeric_limits<double>::quiet_NaN();

    GlobalDimMatrixType const permeability =
        MaterialPropertyLib::formEigenTensor<GlobalDim>(
            medium[MaterialPropertyLib::PropertyType::permeability].value(
                vars, pos, t, dt));

    assert(permeability.rows() == GlobalDim || permeability.rows() == 1);

    bool const is_scalar_permeability = (permeability.size() == 1);

    auto velocity_cache_vectors = MathLib::createZeroedMatrix<
        Eigen::Matrix<double, GlobalDim, Eigen::Dynamic, Eigen::RowMajor>>(
        velocity_cache, GlobalDim, n_integration_points);

    computeDarcyVelocity(is_scalar_permeability, t, dt, local_x, pos,
                         velocity_cache_vectors);

    return velocity_cache;
}

template <typename ShapeFunction, int GlobalDim>
template <typename LaplacianGravityVelocityCalculator,
          typename VelocityCacheType>
void LiquidFlowLocalAssembler<ShapeFunction, GlobalDim>::
    computeProjectedDarcyVelocity(
        const double t, const double dt, std::vector<double> const& local_x,
        ParameterLib::SpatialPosition const& pos,
        VelocityCacheType& darcy_velocity_at_ips) const
{
    auto const local_matrix_size = local_x.size();
    assert(local_matrix_size == ShapeFunction::NPOINTS);

    const auto local_p_vec =
        MathLib::toVector<NodalVectorType>(local_x, local_matrix_size);

    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    auto const& medium = *_process_data.media_map.getMedium(_element.getID());
    auto const& fluid_phase = fluidPhase(medium);

    MaterialPropertyLib::VariableArray vars;
    auto& phase_pressure = medium.hasPhase("Gas") ? vars.gas_phase_pressure
                                                  : vars.liquid_phase_pressure;

    vars.temperature =
        medium[MaterialPropertyLib::PropertyType::reference_temperature]
            .template value<double>(vars, pos, t, dt);

    GlobalDimVectorType const projected_body_force_vector =
        _process_data.element_rotation_matrices[_element.getID()] *
        _process_data.element_rotation_matrices[_element.getID()].transpose() *
        _process_data.specific_body_force;

    auto const& Ns = _process_data.shape_matrix_cache
                         .NsHigherOrder<typename ShapeFunction::MeshElement>();

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        auto const& ip_data = _ip_data[ip];
        auto const& N = Ns[ip];

        phase_pressure = N.dot(local_p_vec);

        auto const [fluid_density, viscosity] =
            MaterialPropertyLib::getFluidDensityAndViscosity(t, dt, pos,
                                                             fluid_phase, vars);

        GlobalDimMatrixType const permeability =
            MaterialPropertyLib::formEigenTensor<GlobalDim>(
                medium[MaterialPropertyLib::PropertyType::permeability].value(
                    vars, pos, t, dt));

        darcy_velocity_at_ips.col(ip) =
            LaplacianGravityVelocityCalculator::calculateVelocity(
                local_p_vec, ip_data, permeability, viscosity, fluid_density,
                projected_body_force_vector, _process_data.has_gravity);
    }
}

template <typename ShapeFunction, int GlobalDim>
void LiquidFlowLocalAssembler<ShapeFunction, GlobalDim>::IsotropicCalculator::
    calculateLaplacianAndGravityTerm(
        Eigen::Map<NodalMatrixType>& local_K,
        Eigen::Map<NodalVectorType>& local_b,
        IntegrationPointData<GlobalDimNodalMatrixType> const& ip_data,
        GlobalDimMatrixType const& permeability_with_density_factor,
        double const mu, double const rho_L,
        GlobalDimVectorType const& specific_body_force, bool const has_gravity)
{
    const double K = permeability_with_density_factor(0, 0) / mu;
    const double fac = K * ip_data.integration_weight;
    local_K.noalias() += fac * ip_data.dNdx.transpose() * ip_data.dNdx;

    if (has_gravity)
    {
        local_b.noalias() +=
            (fac * rho_L) * ip_data.dNdx.transpose() * specific_body_force;
    }
}

template <typename ShapeFunction, int GlobalDim>
Eigen::Matrix<double, GlobalDim, 1>
LiquidFlowLocalAssembler<ShapeFunction, GlobalDim>::IsotropicCalculator::
    calculateVelocity(
        Eigen::Map<const NodalVectorType> const& local_p,
        IntegrationPointData<GlobalDimNodalMatrixType> const& ip_data,
        GlobalDimMatrixType const& permeability, double const mu,
        double const rho_L, GlobalDimVectorType const& specific_body_force,
        bool const has_gravity)
{
    const double K = permeability(0, 0) / mu;
    // Compute the velocity
    GlobalDimVectorType velocity = -K * ip_data.dNdx * local_p;
    // gravity term
    if (has_gravity)
    {
        velocity += (K * rho_L) * specific_body_force;
    }
    return velocity;
}

template <typename ShapeFunction, int GlobalDim>
void LiquidFlowLocalAssembler<ShapeFunction, GlobalDim>::AnisotropicCalculator::
    calculateLaplacianAndGravityTerm(
        Eigen::Map<NodalMatrixType>& local_K,
        Eigen::Map<NodalVectorType>& local_b,
        IntegrationPointData<GlobalDimNodalMatrixType> const& ip_data,
        GlobalDimMatrixType const& permeability_with_density_factor,
        double const mu, double const rho_L,
        GlobalDimVectorType const& specific_body_force, bool const has_gravity)
{
    const double fac = ip_data.integration_weight / mu;
    local_K.noalias() += fac * ip_data.dNdx.transpose() *
                         permeability_with_density_factor * ip_data.dNdx;

    if (has_gravity)
    {
        local_b.noalias() += (fac * rho_L) * ip_data.dNdx.transpose() *
                             permeability_with_density_factor *
                             specific_body_force;
    }
}

template <typename ShapeFunction, int GlobalDim>
Eigen::Matrix<double, GlobalDim, 1>
LiquidFlowLocalAssembler<ShapeFunction, GlobalDim>::AnisotropicCalculator::
    calculateVelocity(
        Eigen::Map<const NodalVectorType> const& local_p,
        IntegrationPointData<GlobalDimNodalMatrixType> const& ip_data,
        GlobalDimMatrixType const& permeability, double const mu,
        double const rho_L, GlobalDimVectorType const& specific_body_force,
        bool const has_gravity)
{
    // Compute the velocity
    GlobalDimVectorType velocity = -permeability * ip_data.dNdx * local_p / mu;

    // gravity term
    if (has_gravity)
    {
        velocity += (rho_L / mu) * permeability * specific_body_force;
    }
    return velocity;
}

}  // namespace LiquidFlow
}  // namespace ProcessLib
