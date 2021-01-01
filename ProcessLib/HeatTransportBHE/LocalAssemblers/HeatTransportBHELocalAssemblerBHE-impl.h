/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "HeatTransportBHELocalAssemblerBHE.h"
#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "NumLib/Fem/InitShapeMatrices.h"

namespace ProcessLib
{
namespace HeatTransportBHE
{
template <typename ShapeFunction, typename IntegrationMethod, typename BHEType>
HeatTransportBHELocalAssemblerBHE<ShapeFunction, IntegrationMethod, BHEType>::
    HeatTransportBHELocalAssemblerBHE(MeshLib::Element const& e,
                                      BHEType const& bhe,
                                      bool const is_axially_symmetric,
                                      unsigned const integration_order,
                                      HeatTransportBHEProcessData& process_data)
    : _process_data(process_data),
      _integration_method(integration_order),
      _bhe(bhe),
      _element_id(e.getID())
{
    // need to make sure that the BHE elements are one-dimensional
    assert(e.getDimension() == 1);

    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    _ip_data.reserve(n_integration_points);
    _secondary_data.N.resize(n_integration_points);

    auto const shape_matrices =
        NumLib::initShapeMatrices<ShapeFunction, ShapeMatricesType,
                                  3 /* GlobalDim */>(e, is_axially_symmetric,
                                                     _integration_method);

    // ip data initialization
    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        auto const& sm = shape_matrices[ip];
        // create the class IntegrationPointDataBHE in place
        _ip_data.push_back(
            {sm.N, sm.dNdx,
             _integration_method.getWeightedPoint(ip).getWeight() *
                 sm.integralMeasure * sm.detJ});

        _secondary_data.N[ip] = sm.N;
    }

    // calculate the element direction vector
    auto const p0 =
        Eigen::Map<Eigen::Vector3d const>(e.getNode(0)->getCoords(), 3);
    auto const p1 =
        Eigen::Map<Eigen::Vector3d const>(e.getNode(1)->getCoords(), 3);

    _element_direction = (p1 - p0).normalized();

    _R_matrix.setZero(bhe_unknowns_size, bhe_unknowns_size);
    _R_pi_s_matrix.setZero(bhe_unknowns_size, soil_temperature_size);
    _R_s_matrix.setZero(soil_temperature_size, soil_temperature_size);
    static constexpr int max_num_thermal_exchange_terms = 5;
    // formulate the local BHE R matrix
    for (int idx_bhe_unknowns = 0; idx_bhe_unknowns < bhe_unknowns;
         idx_bhe_unknowns++)
    {
        typename ShapeMatricesType::template MatrixType<
            single_bhe_unknowns_size, single_bhe_unknowns_size>
            matBHE_loc_R = ShapeMatricesType::template MatrixType<
                single_bhe_unknowns_size,
                single_bhe_unknowns_size>::Zero(single_bhe_unknowns_size,
                                                single_bhe_unknowns_size);
        // Loop over Gauss points
        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            auto const& N = _ip_data[ip].N;
            auto const& w = _ip_data[ip].integration_weight;

            auto const& R = _bhe.thermalResistance(idx_bhe_unknowns);
            // calculate mass matrix for current unknown
            matBHE_loc_R += N.transpose() * N * (1 / R) * w;
        }  // end of loop over integration point

        // The following assembly action is according to Diersch (2013) FEFLOW
        // book please refer to M.127 and M.128 on page 955 and 956
        // The if check is absolutely necessary because
        // (i) In the CXA and CXC case, there are 3 exchange terms,
        // and it is the same as the number of unknowns;
        // (ii) In the 1U case, there are 4 exchange terms,
        // and it is again same as the number of unknowns;
        // (iii) In the 2U case, there are 5 exchange terms,
        // and it is less than the number of unknowns (8).
        if (idx_bhe_unknowns < max_num_thermal_exchange_terms)
        {
            _bhe.template assembleRMatrices<ShapeFunction::NPOINTS>(
                idx_bhe_unknowns, matBHE_loc_R, _R_matrix, _R_pi_s_matrix,
                _R_s_matrix);
        }
    }  // end of loop over BHE unknowns

    // debugging
    // std::string sep =
    //     "\n----------------------------------------\n";
    // Eigen::IOFormat CleanFmt(4, 0, ", ", "\n", "[", "]");
    // std::cout << "_R_matrix: \n" << sep;
    // std::cout << _R_matrix.format(CleanFmt) << sep;
    // std::cout << "_R_s_matrix: \n" << sep;
    // std::cout << _R_s_matrix.format(CleanFmt) << sep;
    // std::cout << "_R_pi_s_matrix: \n" << sep;
    // std::cout << _R_pi_s_matrix.format(CleanFmt) << sep;
}

template <typename ShapeFunction, typename IntegrationMethod, typename BHEType>
void HeatTransportBHELocalAssemblerBHE<ShapeFunction, IntegrationMethod,
                                       BHEType>::
    assemble(
        double const /*t*/, double const /*dt*/,
        std::vector<double> const& /*local_x*/,
        std::vector<double> const& /*local_xdot*/,
        std::vector<double>& local_M_data, std::vector<double>& local_K_data,
        std::vector<double>& /*local_b_data*/)  // local b vector is not touched
{
    auto local_M = MathLib::createZeroedMatrix<BheLocalMatrixType>(
        local_M_data, local_matrix_size, local_matrix_size);
    auto local_K = MathLib::createZeroedMatrix<BheLocalMatrixType>(
        local_K_data, local_matrix_size, local_matrix_size);

    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    auto const& pipe_heat_capacities = _bhe.pipeHeatCapacities();
    auto const& pipe_heat_conductions = _bhe.pipeHeatConductions();
    auto const& pipe_advection_vectors =
        _bhe.pipeAdvectionVectors(_element_direction);
    auto const& cross_section_areas = _bhe.crossSectionAreas();

    // the mass and conductance matrix terms
    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        auto& ip_data = _ip_data[ip];

        auto const& w = ip_data.integration_weight;
        auto const& N = ip_data.N;
        auto const& dNdx = ip_data.dNdx;

        // looping over all unknowns.
        for (int idx_bhe_unknowns = 0; idx_bhe_unknowns < bhe_unknowns;
             idx_bhe_unknowns++)
        {
            // get coefficient of mass from corresponding BHE.
            auto const& mass_coeff = pipe_heat_capacities[idx_bhe_unknowns];
            auto const& lambda = pipe_heat_conductions[idx_bhe_unknowns];
            auto const& advection_vector =
                pipe_advection_vectors[idx_bhe_unknowns];
            auto const& A = cross_section_areas[idx_bhe_unknowns];

            int const single_bhe_unknowns_index =
                bhe_unknowns_index +
                single_bhe_unknowns_size * idx_bhe_unknowns;
            // local M
            local_M
                .template block<single_bhe_unknowns_size,
                                single_bhe_unknowns_size>(
                    single_bhe_unknowns_index, single_bhe_unknowns_index)
                .noalias() += N.transpose() * N * mass_coeff * A * w;

            // local K
            // laplace part
            local_K
                .template block<single_bhe_unknowns_size,
                                single_bhe_unknowns_size>(
                    single_bhe_unknowns_index, single_bhe_unknowns_index)
                .noalias() += dNdx.transpose() * dNdx * lambda * A * w;
            // advection part
            local_K
                .template block<single_bhe_unknowns_size,
                                single_bhe_unknowns_size>(
                    single_bhe_unknowns_index, single_bhe_unknowns_index)
                .noalias() +=
                N.transpose() * advection_vector.transpose() * dNdx * A * w;
        }
    }

    // add the R matrix to local_K
    local_K.template block<bhe_unknowns_size, bhe_unknowns_size>(
        bhe_unknowns_index, bhe_unknowns_index) += _R_matrix;

    // add the R_pi_s matrix to local_K
    local_K
        .template block<bhe_unknowns_size, soil_temperature_size>(
            bhe_unknowns_index, soil_temperature_index)
        .noalias() += _R_pi_s_matrix;
    local_K
        .template block<soil_temperature_size, bhe_unknowns_size>(
            soil_temperature_index, bhe_unknowns_index)
        .noalias() += _R_pi_s_matrix.transpose();

    // add the R_s matrix to local_K
    local_K
        .template block<soil_temperature_size, soil_temperature_size>(
            soil_temperature_index, soil_temperature_index)
        .noalias() += _bhe.number_of_grout_zones * _R_s_matrix;

    // debugging
    // std::string sep = "\n----------------------------------------\n";
    // Eigen::IOFormat CleanFmt(4, 0, ", ", "\n", "[", "]");
    // std::cout << local_K.format(CleanFmt) << sep;
    // std::cout << local_M.format(CleanFmt) << sep;
}
}  // namespace HeatTransportBHE
}  // namespace ProcessLib
