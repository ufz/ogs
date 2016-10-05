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

#ifndef OGS_LIQUIDFLOWLOCALASSEMBLER_IMPL_H
#define OGS_LIQUIDFLOWLOCALASSEMBLER_IMPL_H

#include "LiquidFlowLocalAssembler.h"

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
    _material_properties.setMaterialID(pos);

    const Eigen::MatrixXd& perm =
        _material_properties.getPermeability(t, pos, _element.getDimension());
    // Note: For Inclined 1D in 2D/3D or 2D element in 3D, the first item in
    //  the assert must be changed to perm.rows() == _element->getDimension()
    assert(perm.rows() == GlobalDim || perm.rows() == 1);

    if (perm.size() == 1)  // isotropic or 1D problem.
        local_assemble<IsotropicLaplacianAndGravityTermCalculator>(
            t, local_x, local_M_data, local_K_data, local_b_data, pos, perm);
    else
        local_assemble<AnisotropicLaplacianAndGravityTermCalculator>(
            t, local_x, local_M_data, local_K_data, local_b_data, pos, perm);
}

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
template <typename LaplacianAndGravityTermCalculator>
void LiquidFlowLocalAssembler<ShapeFunction, IntegrationMethod, GlobalDim>::
    local_assemble(double const t, std::vector<double> const& local_x,
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
    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        auto const& sm = _shape_matrices[ip];
        auto const& wp = _integration_method.getWeightedPoint(ip);

        double p = 0.;
        NumLib::shapeFunctionInterpolate(local_x, sm.N, p);
        // TODO : compute _temperature from the heat transport pcs

        const double integration_factor =
            sm.integralMeasure * sm.detJ * wp.getWeight();

        // Assemble mass matrix, M
        local_M.noalias() +=
            _material_properties.getMassCoefficient(
                t, pos, porosity_variable, storage_variable, p, _temperature) *
            sm.N.transpose() * sm.N * integration_factor;

        // Compute density:
        const double rho_g =
            _material_properties.getLiquidDensity(p, _temperature) *
            _gravitational_acceleration;
        // Compute viscosity:
        const double mu = _material_properties.getViscosity(p, _temperature);

        // Assemble Laplacian, K, and RHS by the gravitational term
        LaplacianAndGravityTermCalculator::calculate(
            local_K, local_b, sm, perm, integration_factor, mu, rho_g,
            _gravitational_axis_id);
    }
}

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
void LiquidFlowLocalAssembler<ShapeFunction, IntegrationMethod, GlobalDim>::
    IsotropicLaplacianAndGravityTermCalculator::calculate(
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
    AnisotropicLaplacianAndGravityTermCalculator::calculate(
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
    computeSecondaryVariableConcrete(std::vector<double> const& local_x)
{
    auto const local_matrix_size = local_x.size();
    assert(local_matrix_size == ShapeFunction::NPOINTS);

    const auto local_p_vec =
        MathLib::toVector<NodalVectorType>(local_x, local_matrix_size);

    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    SpatialPosition pos;
    pos.setElementID(_element.getID());
    _material_properties.setMaterialID(pos);
    const Eigen::MatrixXd& perm =
        _material_properties.getPermeability(t, pos, _element.getDimension());

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        auto const& sm = _shape_matrices[ip];
        double p = 0.;
        NumLib::shapeFunctionInterpolate(local_x, sm.N, p);
        // TODO : compute _temperature from the heat transport pcs

        const double rho_g =
            _material_properties.getLiquidDensity(p, _temperature) *
            _gravitational_acceleration;
        // Compute viscosity:
        const double mu = _material_properties.getViscosity(p, _temperature);

        // Assemble Laplacian, K, and RHS by the gravitational term
        if (perm.size() == 1)  // Save time for isotropic permeability.
        {
            //  Use scalar number for isotropic permeability
            //  to save the computation time.
            const double K = perm(0, 0) / mu;
            // Compute the velocity
            GlobalDimVectorType darcy_velocity = -K * sm.dNdx * local_p_vec;
            // gravity term
            if (_gravitational_axis_id >= 0)
                darcy_velocity[GlobalDim - 1] -= K * rho_g;
            for (unsigned d = 0; d < GlobalDim; ++d)
            {
                _darcy_velocities[d][ip] = darcy_velocity[d];
            }
        }
        else
        {
            // Compute the velocity
            GlobalDimVectorType darcy_velocity = -perm * sm.dNdx * local_p_vec / mu;
            if (_gravitational_axis_id >= 0)
            {
                darcy_velocity.noalias() -=  rho_g * perm.col(GlobalDim - 1) / mu;
            }
            for (unsigned d = 0; d < GlobalDim; ++d)
            {
                _darcy_velocities[d][ip] = darcy_velocity[d];
            }
        }
    }

}

}  // end of namespace
}  // end of namespace

#endif /* OGS_LIQUIDFLOWLOCALASSEMBLER_IMPL_H */
