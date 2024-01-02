/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "MaterialLib/MPL/Utils/FormEigenTensor.h"
#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"
#include "NumLib/Function/Interpolation.h"
#include "TwoPhaseFlowWithPrhoLocalAssembler.h"
#include "TwoPhaseFlowWithPrhoProcessData.h"

namespace ProcessLib
{
namespace TwoPhaseFlowWithPrho
{
namespace MPL = MaterialPropertyLib;

template <typename ShapeFunction, int GlobalDim>
void TwoPhaseFlowWithPrhoLocalAssembler<ShapeFunction, GlobalDim>::assemble(
    double const t, double const dt, std::vector<double> const& local_x,
    std::vector<double> const& /*local_x_prev*/,
    std::vector<double>& local_M_data, std::vector<double>& local_K_data,
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

    auto const& medium = *_process_data.media_map.getMedium(_element.getID());
    auto const& liquid_phase = medium.phase("AqueousLiquid");
    auto const& gas_phase = medium.phase("Gas");

    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    ParameterLib::SpatialPosition pos;
    pos.setElementID(_element.getID());
    const int material_id =
        _process_data._material->getMaterialID(pos.getElementID().value());

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        auto const& sm = _shape_matrices[ip];

        double pl_int_pt = 0.;
        double totalrho_int_pt =
            0.;  // total mass density of the light component
        NumLib::shapeFunctionInterpolate(local_x, sm.N, pl_int_pt,
                                         totalrho_int_pt);

        const double temperature = _process_data._temperature(t, pos)[0];

        MPL::VariableArray variables;

        variables.liquid_phase_pressure = pl_int_pt;
        variables.temperature = temperature;
        variables.molar_mass =
            gas_phase.property(MPL::PropertyType::molar_mass)
                .template value<double>(variables, pos, t, dt);

        auto const permeability = MPL::formEigenTensor<GlobalDim>(
            medium.property(MPL::PropertyType::permeability)
                .value(variables, pos, t, dt));
        auto const rho_gas = gas_phase.property(MPL::PropertyType::density)
                                 .template value<double>(variables, pos, t, dt);
        auto const rho_h2o = liquid_phase.property(MPL::PropertyType::density)
                                 .template value<double>(variables, pos, t, dt);

        double& Sw = _ip_data[ip].sw;
        /// Here only consider one component in gas phase
        double const X_h2_nonwet = 1.0;
        double& rho_h2_wet = _ip_data[ip].rho_m;
        double& dSwdP = _ip_data[ip].dsw_dpg;
        double& dSwdrho = _ip_data[ip].dsw_drho;
        double& drhoh2wet = _ip_data[ip].drhom_dpg;
        double& drhoh2wet_drho = _ip_data[ip].drhom_drho;
        if (!_ip_data[ip].mat_property.computeConstitutiveRelation(
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
        double const pc = _process_data._material->getCapillaryPressure(
            material_id, t, pos, pl_int_pt, temperature, Sw);
        variables.capillary_pressure = pc;
        variables.liquid_saturation = Sw;

        double const rho_wet = rho_h2o + rho_h2_wet;
        _saturation[ip] = Sw;
        variables.gas_phase_pressure = _pressure_nonwetting[ip] =
            pl_int_pt + pc;

        // Assemble M matrix
        // nonwetting
        double dPC_dSw =
            _process_data._material->getCapillaryPressureDerivative(
                material_id, t, pos, pl_int_pt, temperature, Sw);

        double const porosity =
            medium.property(MPL::PropertyType::porosity)
                .template value<double>(variables, pos, t, dt);

        Mgx.noalias() += porosity * _ip_data[ip].massOperator;

        Mlp.noalias() += porosity * rho_h2o * dSwdP * _ip_data[ip].massOperator;

        Mlx.noalias() +=
            porosity * (1 + dSwdrho * rho_h2o) * _ip_data[ip].massOperator;
        double const k_rel_gas =
            _process_data._material->getNonwetRelativePermeability(
                t, pos, _pressure_nonwetting[ip], temperature, Sw);
        double const mu_gas =
            gas_phase.property(MPL::PropertyType::viscosity)
                .template value<double>(variables, pos, t, dt);
        double const lambda_gas = k_rel_gas / mu_gas;
        double const diffusion_coeff_component_h2 =
            _process_data._diffusion_coeff_component_b(t, pos)[0];

        // wet
        double const k_rel_wet =
            _process_data._material->getWetRelativePermeability(
                t, pos, pl_int_pt, temperature, Sw);
        auto const mu_wet = liquid_phase.property(MPL::PropertyType::viscosity)
                                .template value<double>(variables, pos, t, dt);
        double const lambda_wet = k_rel_wet / mu_wet;

        laplace_operator.noalias() = sm.dNdx.transpose() * permeability *
                                     sm.dNdx * _ip_data[ip].integration_weight;

        Kgp.noalias() +=
            (rho_gas * X_h2_nonwet * lambda_gas * (1 + dPC_dSw * dSwdP) +
             rho_h2_wet * lambda_wet) *
                laplace_operator +
            (Sw * porosity * diffusion_coeff_component_h2 *
             (rho_h2o / rho_wet) * drhoh2wet) *
                _ip_data[ip].diffusionOperator;
        Kgx.noalias() +=
            (rho_gas * X_h2_nonwet * lambda_gas * dPC_dSw * dSwdrho) *
                laplace_operator +
            (Sw * porosity * diffusion_coeff_component_h2 *
             (rho_h2o / rho_wet) * drhoh2wet_drho) *
                _ip_data[ip].diffusionOperator;
        Klp.noalias() += (rho_gas * lambda_gas * (1 + dPC_dSw * dSwdP) +
                          rho_wet * lambda_wet) *
                         laplace_operator;

        Klx.noalias() +=
            (rho_gas * lambda_gas * dPC_dSw * dSwdrho) * laplace_operator;

        if (_process_data._has_gravity)
        {
            auto const& b = _process_data._specific_body_force;
            Bg.noalias() += (rho_gas * rho_gas * lambda_gas +
                             rho_h2_wet * rho_wet * lambda_wet) *
                            sm.dNdx.transpose() * permeability * b *
                            _ip_data[ip].integration_weight;
            Bl.noalias() += (rho_wet * lambda_wet * rho_wet +
                             rho_gas * rho_gas * lambda_gas) *
                            sm.dNdx.transpose() * permeability * b *
                            _ip_data[ip].integration_weight;

        }  // end of has gravity
    }
    if (_process_data._has_mass_lumping)
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
