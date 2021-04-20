/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

/**
 * common nomenclature
 * --------------primary variable----------------------
 * pn_int_pt    pressure for nonwetting phase at each integration point
 * pc_int_pt    capillary pressure at each integration point
 * --------------secondary variable--------------------
 * Sw wetting               phase saturation
 * dSw_dpc                  derivative of wetting phase saturation with respect
 * to capillary pressure
 * rho_nonwet               density of nonwetting phase
 * drhononwet_dpn           derivative of nonwetting phase density with respect
 *to nonwetting phase pressure
 * rho_wet                  density of wetting phase
 * k_rel_nonwet             relative permeability of nonwetting phase
 * mu_nonwet                viscosity of nonwetting phase
 * lambda_nonwet            mobility of nonwetting phase
 * k_rel_wet                relative permeability of wetting phase
 * mu_wet                   viscosity of wetting phase
 * lambda_wet               mobility of wetting phase
 */
#pragma once

#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/MPL/Utils/FormEigenTensor.h"
#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"
#include "NumLib/Function/Interpolation.h"
#include "TwoPhaseFlowWithPPProcessData.h"

namespace ProcessLib
{
namespace TwoPhaseFlowWithPP
{
namespace MPL = MaterialPropertyLib;

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
void TwoPhaseFlowWithPPLocalAssembler<
    ShapeFunction, IntegrationMethod,
    GlobalDim>::assemble(double const t, double const dt,
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
    auto Mgpc = local_M.template block<nonwet_pressure_size, cap_pressure_size>(
        nonwet_pressure_matrix_index, cap_pressure_matrix_index);

    auto Mlpc = local_M.template block<cap_pressure_size, cap_pressure_size>(
        cap_pressure_matrix_index, cap_pressure_matrix_index);

    NodalMatrixType laplace_operator =
        NodalMatrixType::Zero(ShapeFunction::NPOINTS, ShapeFunction::NPOINTS);

    auto Kgp =
        local_K.template block<nonwet_pressure_size, nonwet_pressure_size>(
            nonwet_pressure_matrix_index, nonwet_pressure_matrix_index);

    auto Klp = local_K.template block<cap_pressure_size, nonwet_pressure_size>(
        cap_pressure_matrix_index, nonwet_pressure_matrix_index);

    auto Klpc = local_K.template block<cap_pressure_size, cap_pressure_size>(
        cap_pressure_matrix_index, cap_pressure_matrix_index);

    auto Bg = local_b.template segment<nonwet_pressure_size>(
        nonwet_pressure_matrix_index);

    auto Bl =
        local_b.template segment<cap_pressure_size>(cap_pressure_matrix_index);

    auto const& medium = *_process_data.media_map->getMedium(_element.getID());
    auto const& liquid_phase = medium.phase("AqueousLiquid");
    auto const& gas_phase = medium.phase("Gas");

    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    ParameterLib::SpatialPosition pos;
    pos.setElementID(_element.getID());

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        pos.setIntegrationPoint(ip);
        auto const& w = _ip_data[ip].integration_weight;
        auto const& dNdx = _ip_data[ip].dNdx;
        auto const& M = _ip_data[ip].massOperator;

        double pc_int_pt = 0.;
        double pn_int_pt = 0.;
        NumLib::shapeFunctionInterpolate(local_x, _ip_data[ip].N, pn_int_pt,
                                         pc_int_pt);

        _pressure_wet[ip] = pn_int_pt - pc_int_pt;
        MPL::VariableArray variables;

        variables[static_cast<int>(MPL::Variable::phase_pressure)] = pn_int_pt;
        variables[static_cast<int>(MPL::Variable::capillary_pressure)] =
            pc_int_pt;
        variables[static_cast<int>(MPL::Variable::temperature)] =
            _process_data.temperature(t, pos)[0];
        variables[static_cast<int>(MPL::Variable::molar_mass)] =
            gas_phase.property(MPL::PropertyType::molar_mass)
                .template value<double>(variables, pos, t, dt);

        auto const rho_nonwet =
            gas_phase.property(MPL::PropertyType::density)
                .template value<double>(variables, pos, t, dt);

        auto const rho_wet = liquid_phase.property(MPL::PropertyType::density)
                                 .template value<double>(variables, pos, t, dt);
        auto& Sw = _saturation[ip];
        Sw = medium.property(MPL::PropertyType::saturation)
                 .template value<double>(variables, pos, t, dt);

        auto const dSw_dpc =
            medium.property(MPL::PropertyType::saturation)
                .template dValue<double>(
                    variables, MPL::Variable::capillary_pressure, pos, t, dt);

        auto const porosity =
            medium.property(MPL::PropertyType::porosity)
                .template value<double>(variables, pos, t, dt);

        auto const drhononwet_dpn =
            gas_phase.property(MPL::PropertyType::density)
                .template dValue<double>(
                    variables, MPL::Variable::phase_pressure, pos, t, dt);

        auto const k_rel_wet =
            medium.property(MPL::PropertyType::relative_permeability)
                .template value<double>(variables, pos, t, dt);
        auto const k_rel_nonwet =
            medium
                .property(
                    MPL::PropertyType::relative_permeability_nonwetting_phase)
                .template value<double>(variables, pos, t, dt);

        auto const mu_nonwet =
            gas_phase.property(MPL::PropertyType::viscosity)
                .template value<double>(variables, pos, t, dt);

        auto const lambda_nonwet = k_rel_nonwet / mu_nonwet;

        auto const mu_wet = liquid_phase.property(MPL::PropertyType::viscosity)
                                .template value<double>(variables, pos, t, dt);

        auto const lambda_wet = k_rel_wet / mu_wet;

        auto const K = MPL::formEigenTensor<GlobalDim>(
            medium.property(MPL::PropertyType::permeability)
                .template value<double>(variables, pos, t, dt));

        Mgp.noalias() += porosity * (1 - Sw) * drhononwet_dpn * M;
        Mgpc.noalias() += -porosity * rho_nonwet * dSw_dpc * M;
        Mlpc.noalias() += porosity * dSw_dpc * rho_wet * M;

        laplace_operator.noalias() = dNdx.transpose() * K * dNdx * w;

        Kgp.noalias() += rho_nonwet * lambda_nonwet * laplace_operator;
        Klp.noalias() += rho_wet * lambda_wet * laplace_operator;
        Klpc.noalias() += -rho_wet * lambda_wet * laplace_operator;

        if (_process_data.has_gravity)
        {
            auto const& b = _process_data.specific_body_force;

            NodalVectorType gravity_operator = dNdx.transpose() * K * b * w;
            Bg.noalias() +=
                rho_nonwet * rho_nonwet * lambda_nonwet * gravity_operator;
            Bl.noalias() += rho_wet * rho_wet * lambda_wet * gravity_operator;
        }  // end of has gravity
    }
    if (_process_data.has_mass_lumping)
    {
        for (unsigned row = 0; row < Mgpc.cols(); row++)
        {
            for (unsigned column = 0; column < Mgpc.cols(); column++)
            {
                if (row != column)
                {
                    Mgpc(row, row) += Mgpc(row, column);
                    Mgpc(row, column) = 0.0;
                    Mgp(row, row) += Mgp(row, column);
                    Mgp(row, column) = 0.0;
                    Mlpc(row, row) += Mlpc(row, column);
                    Mlpc(row, column) = 0.0;
                }
            }
        }
    }  // end of mass-lumping
}

}  // namespace TwoPhaseFlowWithPP
}  // namespace ProcessLib
