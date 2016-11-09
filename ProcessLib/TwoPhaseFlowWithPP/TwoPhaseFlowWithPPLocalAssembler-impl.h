/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef OGS_TWOPHASEFLOWWITHPPLOCALASSEMBLER_IMPL_H
#define OGS_TWOPHASEFLOWWITHPPLOCALASSEMBLER_IMPL_H

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

    NodalMatrixType mass_operator;
    mass_operator.setZero(ShapeFunction::NPOINTS, ShapeFunction::NPOINTS);

    auto Mgp =
        local_M.template block<nonwet_pressure_size, nonwet_pressure_size>(
            nonwet_pressure_matrix_index, nonwet_pressure_matrix_index);
    auto Mgpc = local_M.template block<nonwet_pressure_size, cap_pressure_size>(
        nonwet_pressure_matrix_index, cap_pressure_matrix_index);

    auto Mlp = local_M.template block<cap_pressure_size, nonwet_pressure_size>(
        cap_pressure_matrix_index, nonwet_pressure_matrix_index);

    auto Mlpc = local_M.template block<cap_pressure_size, cap_pressure_size>(
        cap_pressure_matrix_index, cap_pressure_matrix_index);

    NodalMatrixType laplace_operator;
    laplace_operator.setZero(ShapeFunction::NPOINTS, ShapeFunction::NPOINTS);

    auto Kgp =
        local_K.template block<nonwet_pressure_size, nonwet_pressure_size>(
            nonwet_pressure_matrix_index, nonwet_pressure_matrix_index);

    auto Kgpc = local_K.template block<nonwet_pressure_size, cap_pressure_size>(
        nonwet_pressure_matrix_index, cap_pressure_matrix_index);

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
    _process_data._material->setMaterialID(pos);

    const Eigen::MatrixXd& perm = _process_data._material->getPermeability(
        t, pos, _element.getDimension());
    assert(perm.rows() == GlobalDim || perm.rows() == 1);
    GlobalDimMatrixType permeability =
        GlobalDimMatrixType::Zero(GlobalDim, GlobalDim);
    if (perm.rows() == GlobalDim)
        permeability = perm;
    else if (perm.rows() == 1)
        permeability.diagonal().setConstant(perm(0, 0));
    MathLib::PiecewiseLinearInterpolation const& interpolated_Pc =
        _process_data._interpolated_Pc;
    MathLib::PiecewiseLinearInterpolation const& interpolated_Kr_wet =
        _process_data._interpolated_Kr_wet;
    MathLib::PiecewiseLinearInterpolation const& interpolated_Kr_nonwet =
        _process_data._interpolated_Kr_nonwet;

    // Note: currently only isothermal case is considered, so the temperature is
    // assumed to be const
    // the variation of temperature will be taken into account in future
    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        auto const& sm = _shape_matrices[ip];

        double pc_int_pt = 0.;
        double pg_int_pt = 0.;
        NumLib::shapeFunctionInterpolate(local_x, sm.N, pg_int_pt, pc_int_pt);

        _pressure_wetting[ip] = pg_int_pt - pc_int_pt;

        auto const& wp = _integration_method.getWeightedPoint(ip);
        const double integration_factor =
            sm.integralMeasure * sm.detJ * wp.getWeight();

        double const rho_gas =
            _process_data._material->getGasDensity(pg_int_pt, _temperature);
        double const rho_w = _process_data._material->getLiquidDensity(
            _pressure_wetting[ip], _temperature);

        double const Sw =
            (pc_int_pt < 0) ? 1.0 : interpolated_Pc.getValue(pc_int_pt);

        _saturation[ip] = Sw;
        double dSwdPc = interpolated_Pc.getDerivative(pc_int_pt);
        if (pc_int_pt > interpolated_Pc.getSupportMax())
            dSwdPc =
                interpolated_Pc.getDerivative(interpolated_Pc.getSupportMax());
        else if (pc_int_pt < interpolated_Pc.getSupportMin())
            dSwdPc =
                interpolated_Pc.getDerivative(interpolated_Pc.getSupportMin());

        double const porosity = _process_data._material->getPorosity(
            t, pos, pg_int_pt, _temperature, 0);

        // Assemble M matrix
        // nonwetting
        double const drhogas_dpg = _process_data._material->getDerivGasDensity(
            pg_int_pt, _temperature);

        mass_operator.noalias() = sm.N.transpose() * sm.N * integration_factor;

        Mgp.noalias() += porosity * (1 - Sw) * drhogas_dpg * mass_operator;
        Mgpc.noalias() += -porosity * rho_gas * dSwdPc * mass_operator;
        Mlp.noalias() += 0.0 * mass_operator;
        Mlpc.noalias() += porosity * dSwdPc * rho_w * mass_operator;

        // Assemble M matrix
        // nonwet
        double const k_rel_G = interpolated_Kr_nonwet.getValue(Sw);
        double const mu_gas =
            _process_data._material->getGasViscosity(pg_int_pt, _temperature);
        double const lambda_G = k_rel_G / mu_gas;

        // wet
        double const k_rel_L = interpolated_Kr_wet.getValue(Sw);
        double const mu_liquid = _process_data._material->getLiquidViscosity(
            _pressure_wetting[ip], _temperature);
        double const lambda_L = k_rel_L / mu_liquid;

        laplace_operator.noalias() =
            sm.dNdx.transpose() * permeability * sm.dNdx * integration_factor;

        Kgp.noalias() += rho_gas * lambda_G * laplace_operator;
        Kgpc.noalias() += 0.0 * laplace_operator;
        Klp.noalias() += rho_w * lambda_L * laplace_operator;
        Klpc.noalias() += -rho_w * lambda_L * laplace_operator;

        if (_process_data._has_gravity)
        {
            auto const& b = _process_data._specific_body_force;
            Bg.noalias() += rho_gas * rho_gas * lambda_G * sm.dNdx.transpose() *
                            permeability * b * integration_factor;
            Bl.noalias() += rho_w * rho_w * lambda_L * sm.dNdx.transpose() *
                            permeability * b * integration_factor;

        }  // end of has gravity
    }      // end of GP
    if (_process_data._has_mass_lumping)
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
                    Mlp(row, row) += Mlp(row, column);
                    Mlp(row, column) = 0.0;
                }
            }
        }
    }  // end of mass-lumping
}

}  // end of namespace
}  // end of namespace

#endif
