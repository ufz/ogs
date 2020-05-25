/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "TwoPhaseFlowWithPrhoLocalAssembler.h"

#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"
#include "NumLib/Function/Interpolation.h"
#include "TwoPhaseFlowWithPrhoProcessData.h"

namespace ProcessLib
{
namespace TwoPhaseFlowWithPrho
{
template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
void TwoPhaseFlowWithPrhoLocalAssembler<
    ShapeFunction, IntegrationMethod,
    GlobalDim>::assemble(double const t, double const /*dt*/,
                         std::vector<double> const& local_x,
                         std::vector<double> const& /*local_xdot*/,
                         std::vector<double>& local_M_data,
                         std::vector<double>& local_K_data,
                         std::vector<double>& local_b_data)
{
    auto const local_matrix_size = local_x.size();

    assert(local_matrix_size == ShapeFunction::NPOINTS * NUM_NODAL_DOF);

    auto local_M = MathLib::createZeroedMatrix<LocalMatrixType>(
        local_M_data, local_matrix_size, local_matrix_size);
    auto local_K = MathLib::createZeroedMatrix<LocalMatrixType>(
        local_K_data, local_matrix_size, local_matrix_size);
    auto local_b = MathLib::createZeroedVector<LocalVectorType>(
        local_b_data, local_matrix_size);

    auto Mgp =
        local_M.template block<nonwet_pressure_size, nonwet_pressure_size>(
            nonwet_pressure_matrix_index, nonwet_pressure_matrix_index);
    auto Mgx = local_M.template block<nonwet_pressure_size, cap_pressure_size>(
        nonwet_pressure_matrix_index, cap_pressure_matrix_index);

    auto Mlp = local_M.template block<cap_pressure_size, nonwet_pressure_size>(
        cap_pressure_matrix_index, nonwet_pressure_matrix_index);

    auto Mlx = local_M.template block<cap_pressure_size, cap_pressure_size>(
        cap_pressure_matrix_index, cap_pressure_matrix_index);

    NodalMatrixType laplace_operator =
        NodalMatrixType::Zero(ShapeFunction::NPOINTS, ShapeFunction::NPOINTS);

    auto Kgp =
        local_K.template block<nonwet_pressure_size, nonwet_pressure_size>(
            nonwet_pressure_matrix_index, nonwet_pressure_matrix_index);

    auto Kgx = local_K.template block<nonwet_pressure_size, cap_pressure_size>(
        nonwet_pressure_matrix_index, cap_pressure_matrix_index);

    auto Klp = local_K.template block<cap_pressure_size, nonwet_pressure_size>(
        cap_pressure_matrix_index, nonwet_pressure_matrix_index);

    auto Klx = local_K.template block<cap_pressure_size, cap_pressure_size>(
        cap_pressure_matrix_index, cap_pressure_matrix_index);

    auto Bg = local_b.template segment<nonwet_pressure_size>(
        nonwet_pressure_matrix_index);

    auto Bl =
        local_b.template segment<cap_pressure_size>(cap_pressure_matrix_index);

    unsigned const n_integration_points =
        integration_method_.getNumberOfPoints();

    ParameterLib::SpatialPosition pos;
    pos.setElementID(element_.getID());
    const int material_id =
        process_data_.material_->getMaterialID(pos.getElementID().get());

    const Eigen::MatrixXd& perm = process_data_.material_->getPermeability(
        material_id, t, pos, element_.getDimension());
    assert(perm.rows() == element_.getDimension() || perm.rows() == 1);
    GlobalDimMatrixType permeability = GlobalDimMatrixType::Zero(
        element_.getDimension(), element_.getDimension());
    if (perm.rows() == element_.getDimension())
    {
        permeability = perm;
    }
    else if (perm.rows() == 1)
    {
        permeability.diagonal().setConstant(perm(0, 0));
    }

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        auto const& sm = shape_matrices_[ip];

        double pl_int_pt = 0.;
        double totalrho_int_pt =
            0.;  // total mass density of the light component
        NumLib::shapeFunctionInterpolate(local_x, sm.N, pl_int_pt,
                                         totalrho_int_pt);

        const double temperature = process_data_.temperature_(t, pos)[0];

        double const rho_gas =
            process_data_.material_->getGasDensity(pl_int_pt, temperature);
        double const rho_h2o =
            process_data_.material_->getLiquidDensity(pl_int_pt, temperature);

        double& Sw = ip_data_[ip].sw;
        /// Here only consider one component in gas phase
        double const X_h2_nonwet = 1.0;
        double& rho_h2_wet = ip_data_[ip].rho_m;
        double& dSwdP = ip_data_[ip].dsw_dpg;
        double& dSwdrho = ip_data_[ip].dsw_drho;
        double& drhoh2wet = ip_data_[ip].drhom_dpg;
        double& drhoh2wet_drho = ip_data_[ip].drhom_drho;
        if (!ip_data_[ip].mat_property.computeConstitutiveRelation(
                t,
                pos,
                material_id,
                pl_int_pt,
                totalrho_int_pt,
                temperature,
                Sw,
                rho_h2_wet,
                dSwdP,
                dSwdrho,
                drhoh2wet,
                drhoh2wet_drho))
        {
            OGS_FATAL("Computation of local constitutive relation failed.");
        }
        double const pc = process_data_.material_->getCapillaryPressure(
            material_id, t, pos, pl_int_pt, temperature, Sw);

        double const rho_wet = rho_h2o + rho_h2_wet;
        saturation_[ip] = Sw;
        pressure_nonwetting_[ip] = pl_int_pt + pc;

        // Assemble M matrix
        // nonwetting
        double dPC_dSw =
            process_data_.material_->getCapillaryPressureDerivative(
                material_id, t, pos, pl_int_pt, temperature, Sw);

        double const porosity = process_data_.material_->getPorosity(
            material_id, t, pos, pl_int_pt, temperature, 0);

        Mgx.noalias() += porosity * ip_data_[ip].massOperator;

        Mlp.noalias() += porosity * rho_h2o * dSwdP * ip_data_[ip].massOperator;

        Mlx.noalias() +=
            porosity * (1 + dSwdrho * rho_h2o) * ip_data_[ip].massOperator;
        double const k_rel_gas =
            process_data_.material_->getNonwetRelativePermeability(
                t, pos, pressure_nonwetting_[ip], temperature, Sw);
        double const mu_gas = process_data_.material_->getGasViscosity(
            pressure_nonwetting_[ip], temperature);
        double const lambda_gas = k_rel_gas / mu_gas;
        double const diffusion_coeff_component_h2 =
            process_data_.diffusion_coeff_component_b_(t, pos)[0];

        // wet
        double const k_rel_wet =
            process_data_.material_->getWetRelativePermeability(
                t, pos, pl_int_pt, temperature, Sw);
        double const mu_wet =
            process_data_.material_->getLiquidViscosity(pl_int_pt, temperature);
        double const lambda_wet = k_rel_wet / mu_wet;

        laplace_operator.noalias() = sm.dNdx.transpose() * permeability *
                                     sm.dNdx * ip_data_[ip].integration_weight;

        Kgp.noalias() +=
            (rho_gas * X_h2_nonwet * lambda_gas * (1 + dPC_dSw * dSwdP) +
             rho_h2_wet * lambda_wet) *
                laplace_operator +
            (Sw * porosity * diffusion_coeff_component_h2 *
             (rho_h2o / rho_wet) * drhoh2wet) *
                ip_data_[ip].diffusionOperator;
        Kgx.noalias() +=
            (rho_gas * X_h2_nonwet * lambda_gas * dPC_dSw * dSwdrho) *
                laplace_operator +
            (Sw * porosity * diffusion_coeff_component_h2 *
             (rho_h2o / rho_wet) * drhoh2wet_drho) *
                ip_data_[ip].diffusionOperator;
        Klp.noalias() += (rho_gas * lambda_gas * (1 + dPC_dSw * dSwdP) +
                          rho_wet * lambda_wet) *
                         laplace_operator;

        Klx.noalias() +=
            (rho_gas * lambda_gas * dPC_dSw * dSwdrho) * laplace_operator;

        if (process_data_.has_gravity_)
        {
            auto const& b = process_data_.specific_body_force_;
            Bg.noalias() += (rho_gas * rho_gas * lambda_gas +
                             rho_h2_wet * rho_wet * lambda_wet) *
                            sm.dNdx.transpose() * permeability * b *
                            ip_data_[ip].integration_weight;
            Bl.noalias() += (rho_wet * lambda_wet * rho_wet +
                             rho_gas * rho_gas * lambda_gas) *
                            sm.dNdx.transpose() * permeability * b *
                            ip_data_[ip].integration_weight;

        }  // end of has gravity
    }
    if (process_data_.has_mass_lumping_)
    {
        for (unsigned row = 0; row < Mgp.cols(); row++)
        {
            for (unsigned column = 0; column < Mgp.cols(); column++)
            {
                if (row != column)
                {
                    Mgx(row, row) += Mgx(row, column);
                    Mgx(row, column) = 0.0;
                    Mgp(row, row) += Mgp(row, column);
                    Mgp(row, column) = 0.0;
                    Mlx(row, row) += Mlx(row, column);
                    Mlx(row, column) = 0.0;
                    Mlp(row, row) += Mlp(row, column);
                    Mlp(row, column) = 0.0;
                }
            }
        }
    }  // end of mass-lumping
}

}  // namespace TwoPhaseFlowWithPrho
}  // namespace ProcessLib
