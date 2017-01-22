/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
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
* temperature              capillary pressure
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

#include "TwoPhaseFlowWithPPLocalAssembler.h"

#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"
#include "NumLib/Function/Interpolation.h"
#include "TwoPhaseFlowWithPPProcessData.h"

namespace ProcessLib
{
namespace TwoPhaseFlowWithPP
{
template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
void TwoPhaseFlowWithPPLocalAssembler<
    ShapeFunction, IntegrationMethod,
    GlobalDim>::assemble(double const t, std::vector<double> const& local_x,
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

    NodalMatrixType laplace_operator;
    laplace_operator.setZero(ShapeFunction::NPOINTS, ShapeFunction::NPOINTS);

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

    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    SpatialPosition pos;
    pos.setElementID(_element.getID());
    const int material_id =
        _process_data.material->getMaterialID(pos.getElementID().get());

    const Eigen::MatrixXd& perm = _process_data.material->getPermeability(
        material_id, t, pos, _element.getDimension());
    assert(perm.rows() == _element.getDimension() || perm.rows() == 1);
    GlobalDimMatrixType permeability = GlobalDimMatrixType::Zero(
        _element.getDimension(), _element.getDimension());
    if (perm.rows() == _element.getDimension())
        permeability = perm;
    else if (perm.rows() == 1)
        permeability.diagonal().setConstant(perm(0, 0));

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        double pc_int_pt = 0.;
        double pn_int_pt = 0.;
        NumLib::shapeFunctionInterpolate(local_x, _ip_data[ip].N, pn_int_pt,
                                         pc_int_pt);

        _pressure_wet[ip] = pn_int_pt - pc_int_pt;

        const double temperature = _process_data.temperature(t, pos)[0];
        double const rho_nonwet =
            _process_data.material->getGasDensity(pn_int_pt, temperature);
        double const rho_wet = _process_data.material->getLiquidDensity(
            _pressure_wet[ip], temperature);

        double const Sw = _process_data.material->getSaturation(
            material_id, t, pos, pn_int_pt, temperature, pc_int_pt);

        _saturation[ip] = Sw;

        double dSw_dpc = _process_data.material->getSaturationDerivative(
            material_id, t, pos, pn_int_pt, temperature, Sw);

        double const porosity = _process_data.material->getPorosity(
            material_id, t, pos, pn_int_pt, temperature, 0);

        // Assemble M matrix
        // nonwetting
        double const drhononwet_dpn =
            _process_data.material->getGasDensityDerivative(pn_int_pt,
                                                            temperature);

        Mgp.noalias() +=
            porosity * (1 - Sw) * drhononwet_dpn * _ip_data[ip].massOperator;
        Mgpc.noalias() +=
            -porosity * rho_nonwet * dSw_dpc * _ip_data[ip].massOperator;

        Mlpc.noalias() +=
            porosity * dSw_dpc * rho_wet * _ip_data[ip].massOperator;

        // nonwet
        double const k_rel_nonwet =
            _process_data.material->getNonwetRelativePermeability(
                t, pos, pn_int_pt, temperature, Sw);
        double const mu_nonwet =
            _process_data.material->getGasViscosity(pn_int_pt, temperature);
        double const lambda_nonwet = k_rel_nonwet / mu_nonwet;

        // wet
        double const k_rel_wet =
            _process_data.material->getWetRelativePermeability(
                t, pos, _pressure_wet[ip], temperature, Sw);
        double const mu_wet = _process_data.material->getLiquidViscosity(
            _pressure_wet[ip], temperature);
        double const lambda_wet = k_rel_wet / mu_wet;

        laplace_operator.noalias() = _ip_data[ip].dNdx.transpose() *
                                     permeability * _ip_data[ip].dNdx *
                                     _ip_data[ip].integration_weight;

        Kgp.noalias() += rho_nonwet * lambda_nonwet * laplace_operator;

        Klp.noalias() += rho_wet * lambda_wet * laplace_operator;
        Klpc.noalias() += -rho_wet * lambda_wet * laplace_operator;

        if (_process_data.has_gravity)
        {
            auto const& b = _process_data.specific_body_force;

            NodalVectorType gravity_operator = _ip_data[ip].dNdx.transpose() *
                                               permeability * b *
                                               _ip_data[ip].integration_weight;
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

}  // end of namespace
}  // end of namespace
