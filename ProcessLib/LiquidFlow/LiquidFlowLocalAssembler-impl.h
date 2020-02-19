/**
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
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

#include "MaterialLib/MPL/MaterialSpatialDistributionMap.h"
#include "NumLib/Function/Interpolation.h"

namespace ProcessLib
{
namespace LiquidFlow
{
template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
void LiquidFlowLocalAssembler<ShapeFunction, IntegrationMethod, GlobalDim>::
    assemble(double const t, double const dt,
             std::vector<double> const& local_x,
             std::vector<double> const& /*local_xdot*/,
             std::vector<double>& local_M_data,
             std::vector<double>& local_K_data,
             std::vector<double>& local_b_data)
{
    ParameterLib::SpatialPosition pos;
    pos.setElementID(_element.getID());
    const int material_id = _material_properties.getMaterialID(pos);

    double const pressure = std::nan("");
    const Eigen::MatrixXd& permeability = _material_properties.getPermeability(
        material_id, t, pos, _element.getDimension(), pressure,
        _reference_temperature);
    // Note: For Inclined 1D in 2D/3D or 2D element in 3D, the first item in
    //  the assert must be changed to permeability.rows() ==
    //  _element->getDimension()
    assert(permeability.rows() == GlobalDim || permeability.rows() == 1);

    if (permeability.size() == 1)
    {  // isotropic or 1D problem.
        assembleMatrixAndVector<IsotropicCalculator>(
            material_id, t, dt, local_x, local_M_data, local_K_data,
            local_b_data);
    }
    else
    {
        assembleMatrixAndVector<AnisotropicCalculator>(
            material_id, t, dt, local_x, local_M_data, local_K_data,
            local_b_data);
    }
}

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
Eigen::Vector3d
LiquidFlowLocalAssembler<ShapeFunction, IntegrationMethod, GlobalDim>::getFlux(
    MathLib::Point3d const& p_local_coords, double const t,
    std::vector<double> const& local_x) const
{
    // TODO (tf) Temporary value not used by current material models. Need
    // extension of getFlux interface
    double const dt = std::numeric_limits<double>::quiet_NaN();

    // eval dNdx and invJ at p
    auto const fe =
        NumLib::createIsoparametricFiniteElement<ShapeFunction,
                                                 ShapeMatricesType>(_element);

    typename ShapeMatricesType::ShapeMatrices shape_matrices(
        ShapeFunction::DIM, GlobalDim, ShapeFunction::NPOINTS);

    // Note: Axial symmetry is set to false here, because we only need dNdx
    // here, which is not affected by axial symmetry.
    fe.computeShapeFunctions(p_local_coords.getCoords(), shape_matrices,
                             GlobalDim, false);

    // create pos object to access the correct media property
    ParameterLib::SpatialPosition pos;
    pos.setElementID(_element.getID());
    const int material_id = _material_properties.getMaterialID(pos);

    double pressure = 0.0;
    NumLib::shapeFunctionInterpolate(local_x, shape_matrices.N, pressure);
    const Eigen::MatrixXd& permeability = _material_properties.getPermeability(
        material_id, t, pos, _element.getDimension(), pressure,
        _reference_temperature);
    const double mu =
        _material_properties.getViscosity(pressure, _reference_temperature);

    Eigen::Vector3d flux(0.0, 0.0, 0.0);
    flux.head<GlobalDim>() =
        -permeability / mu * shape_matrices.dNdx *
        Eigen::Map<const NodalVectorType>(local_x.data(), local_x.size());

    return flux;
}

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
template <typename LaplacianGravityVelocityCalculator>
void LiquidFlowLocalAssembler<ShapeFunction, IntegrationMethod, GlobalDim>::
    assembleMatrixAndVector(const int material_id, double const t,
                            double const dt, std::vector<double> const& local_x,
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

    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    ParameterLib::SpatialPosition pos;
    pos.setElementID(_element.getID());

    // TODO: The following two variables should be calculated inside the
    //       the integration loop for non-constant porosity and storage models.
    double porosity_variable = 0.;
    double storage_variable = 0.;
    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        auto const& ip_data = _ip_data[ip];

        double p = 0.;
        NumLib::shapeFunctionInterpolate(local_x, ip_data.N, p);

        // Assemble mass matrix, M
        local_M.noalias() += _material_properties.getMassCoefficient(
                                 material_id, t, pos, porosity_variable,
                                 storage_variable, p, _reference_temperature) *
                             ip_data.N.transpose() * ip_data.N *
                             ip_data.integration_weight;

        // Compute density:
        const double rho_g =
            _material_properties.getLiquidDensity(p, _reference_temperature) *
            _gravitational_acceleration;
        // Compute viscosity:
        const double mu =
            _material_properties.getViscosity(p, _reference_temperature);

        pos.setIntegrationPoint(ip);
        auto const& permeability = _material_properties.getPermeability(
            material_id, t, pos, _element.getDimension(), p,
            _reference_temperature);

        // Assemble Laplacian, K, and RHS by the gravitational term
        LaplacianGravityVelocityCalculator::calculateLaplacianAndGravityTerm(
            local_K, local_b, ip_data, permeability, mu, rho_g,
            _gravitational_axis_id);
    }
}

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
std::vector<double> const&
LiquidFlowLocalAssembler<ShapeFunction, IntegrationMethod, GlobalDim>::
    getIntPtDarcyVelocity(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& velocity_cache) const
{
    constexpr int process_id = 0;
    auto const indices =
        NumLib::getIndices(_element.getID(), *dof_table[process_id]);
    auto const local_x = x[process_id]->get(indices);
    auto const num_intpts = _integration_method.getNumberOfPoints();
    velocity_cache.clear();
    auto velocity_cache_vectors = MathLib::createZeroedMatrix<
        Eigen::Matrix<double, GlobalDim, Eigen::Dynamic, Eigen::RowMajor>>(
        velocity_cache, GlobalDim, num_intpts);

    ParameterLib::SpatialPosition pos;
    pos.setElementID(_element.getID());
    const int material_id = _material_properties.getMaterialID(pos);
    // evaluate the permeability to distinguish which computeDarcyVelocity
    // method should be used
    double const pressure = std::nan("");
    const Eigen::MatrixXd& permeability = _material_properties.getPermeability(
        material_id, t, pos, _element.getDimension(), pressure,
        _reference_temperature);

    // Note: For Inclined 1D in 2D/3D or 2D element in 3D, the first item in
    //  the assert must be changed to perm.rows() == _element->getDimension()
    assert(permeability.rows() == GlobalDim || permeability.rows() == 1);

    if (permeability.size() == 1)
    {  // isotropic or 1D problem.
        computeDarcyVelocityLocal<IsotropicCalculator>(
            material_id, t, dt, local_x, pos, velocity_cache_vectors);
    }
    else
    {
        computeDarcyVelocityLocal<AnisotropicCalculator>(
            material_id, t, dt, local_x, pos, velocity_cache_vectors);
    }
    return velocity_cache;
}

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
template <typename LaplacianGravityVelocityCalculator>
void LiquidFlowLocalAssembler<ShapeFunction, IntegrationMethod, GlobalDim>::
    computeDarcyVelocityLocal(
        const int material_id, const double t, const double dt,
        std::vector<double> const& local_x,
        ParameterLib::SpatialPosition const& pos,
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
        auto const& ip_data = _ip_data[ip];
        double p = 0.;
        NumLib::shapeFunctionInterpolate(local_x, ip_data.N, p);

        const double rho_g =
            _material_properties.getLiquidDensity(p, _reference_temperature) *
            _gravitational_acceleration;
        // Compute viscosity:
        const double mu =
            _material_properties.getViscosity(p, _reference_temperature);

        auto const& permeability = _material_properties.getPermeability(
            material_id, t, pos, _element.getDimension(), p,
            _reference_temperature);
        LaplacianGravityVelocityCalculator::calculateVelocity(
            ip, local_p_vec, ip_data, permeability, mu, rho_g,
            _gravitational_axis_id, darcy_velocity_at_ips);
    }
}

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
void LiquidFlowLocalAssembler<ShapeFunction, IntegrationMethod, GlobalDim>::
    IsotropicCalculator::calculateLaplacianAndGravityTerm(
        Eigen::Map<NodalMatrixType>& local_K,
        Eigen::Map<NodalVectorType>& local_b,
        IntegrationPointData<NodalRowVectorType,
                             GlobalDimNodalMatrixType> const& ip_data,
        Eigen::MatrixXd const& permeability, double const mu,
        double const rho_g, int const gravitational_axis_id)
{
    const double K = permeability(0, 0) / mu;
    const double fac = K * ip_data.integration_weight;
    local_K.noalias() += fac * ip_data.dNdx.transpose() * ip_data.dNdx;

    if (gravitational_axis_id >= 0)
    {
        // Note: Since permeability, K, is a scalar number in this function,
        // (gradN)*K*V becomes K*(grad N)*V. With the gravity vector of V,
        // the simplification of (grad N)*V is the gravitational_axis_id th
        // column of the transpose of (grad N) multiplied with rho_g.
        local_b.noalias() -=
            fac * ip_data.dNdx.transpose().col(gravitational_axis_id) * rho_g;
    }
}

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
void LiquidFlowLocalAssembler<ShapeFunction, IntegrationMethod, GlobalDim>::
    IsotropicCalculator::calculateVelocity(
        unsigned const ip, Eigen::Map<const NodalVectorType> const& local_p,
        IntegrationPointData<NodalRowVectorType,
                             GlobalDimNodalMatrixType> const& ip_data,
        Eigen::MatrixXd const& permeability, double const mu,
        double const rho_g, int const gravitational_axis_id,
        MatrixOfVelocityAtIntegrationPoints& darcy_velocity_at_ips)
{
    const double K = permeability(0, 0) / mu;
    // Compute the velocity
    darcy_velocity_at_ips.col(ip).noalias() = -K * ip_data.dNdx * local_p;
    // gravity term
    if (gravitational_axis_id >= 0)
    {
        darcy_velocity_at_ips.col(ip)[gravitational_axis_id] -= K * rho_g;
    }
}

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
void LiquidFlowLocalAssembler<ShapeFunction, IntegrationMethod, GlobalDim>::
    AnisotropicCalculator::calculateLaplacianAndGravityTerm(
        Eigen::Map<NodalMatrixType>& local_K,
        Eigen::Map<NodalVectorType>& local_b,
        IntegrationPointData<NodalRowVectorType,
                             GlobalDimNodalMatrixType> const& ip_data,
        Eigen::MatrixXd const& permeability, double const mu,
        double const rho_g, int const gravitational_axis_id)
{
    const double fac = ip_data.integration_weight / mu;
    local_K.noalias() +=
        fac * ip_data.dNdx.transpose() * permeability * ip_data.dNdx;
    if (gravitational_axis_id >= 0)
    {
        // Note: permeability * gravity_vector = rho_g *
        // permeability.col(gravitational_axis_id)
        //       i.e the result is the gravitational_axis_id th column of
        //       permeability multiplied with rho_g.
        local_b.noalias() -= fac * rho_g * ip_data.dNdx.transpose() *
                             permeability.col(gravitational_axis_id);
    }
}

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
void LiquidFlowLocalAssembler<ShapeFunction, IntegrationMethod, GlobalDim>::
    AnisotropicCalculator::calculateVelocity(
        unsigned const ip, Eigen::Map<const NodalVectorType> const& local_p,
        IntegrationPointData<NodalRowVectorType,
                             GlobalDimNodalMatrixType> const& ip_data,
        Eigen::MatrixXd const& permeability, double const mu,
        double const rho_g, int const gravitational_axis_id,
        MatrixOfVelocityAtIntegrationPoints& darcy_velocity_at_ips)
{
    // Compute the velocity
    darcy_velocity_at_ips.col(ip).noalias() =
        -permeability * ip_data.dNdx * local_p / mu;
    if (gravitational_axis_id >= 0)
    {
        darcy_velocity_at_ips.col(ip).noalias() -=
            rho_g * permeability.col(gravitational_axis_id) / mu;
    }
}

}  // namespace LiquidFlow
}  // namespace ProcessLib
