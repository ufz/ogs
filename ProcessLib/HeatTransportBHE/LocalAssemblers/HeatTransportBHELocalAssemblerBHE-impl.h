// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <algorithm>

#include "BaseLib/Error.h"
#include "HeatTransportBHELocalAssemblerBHE.h"
#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "NumLib/Fem/InitShapeMatrices.h"

namespace ProcessLib
{
namespace HeatTransportBHE
{
template <typename ShapeFunction, typename BHEType>
HeatTransportBHELocalAssemblerBHE<ShapeFunction, BHEType>::
    HeatTransportBHELocalAssemblerBHE(
        MeshLib::Element const& e,
        NumLib::GenericIntegrationMethod const& integration_method,
        BHEType const& bhe,
        bool const is_axially_symmetric,
        HeatTransportBHEProcessData& process_data,
        BHEMeshData const& bhe_mesh_data)
    : _process_data(process_data),
      _integration_method(integration_method),
      _bhe(bhe),
      _element_id(e.getID()),
      _bhe_mesh_data(bhe_mesh_data)
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
    auto const& p0 = e.getNode(0)->asEigenVector3d();
    auto const& p1 = e.getNode(1)->asEigenVector3d();

    _element_direction = (p1 - p0).normalized();

    auto const section_it =
        _bhe_mesh_data.BHE_element_section_indices.find(_element_id);
    if (section_it == _bhe_mesh_data.BHE_element_section_indices.end())
    {
        OGS_FATAL(
            "Could not read BHE element section index for element id "
            "{:d}.",
            _element_id);
    }

    _section_index = section_it->second;
    if (_section_index < 0)
    {
        OGS_FATAL(
            "Invalid BHE section index for element id {:d}. Check BHE mesh "
            "data initialisation.",
            _element_id);
    }

    _R_matrix.setZero(bhe_unknowns_size, bhe_unknowns_size);
    _R_pi_s_matrix.setZero(bhe_unknowns_size, soil_temperature_size);
    _R_s_matrix.setZero(soil_temperature_size, soil_temperature_size);
    static constexpr int max_num_thermal_exchange_terms = 5;
    // Formulate the local BHE R matrix.
    // Only unknowns with thermal exchange terms need resistance assembly.
    // In CXA/CXC there are 3 exchange terms (= number of unknowns),
    // in 1U there are 4 (= number of unknowns),
    // in 2U there are 5 but 8 unknowns — unknowns 5-7 (extra grout zones)
    // have no exchange terms. See Diersch (2013) FEFLOW, M.127-M.128.
    for (int idx_bhe_unknowns = 0;
         idx_bhe_unknowns <
         std::min(bhe_unknowns, max_num_thermal_exchange_terms);
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

            // Get thermal resistance for this element's section
            auto const& R = _bhe.thermalResistanceAtSection(idx_bhe_unknowns,
                                                            _section_index);
            // calculate mass matrix for current unknown
            matBHE_loc_R += N.transpose() * N * (1 / R) * w;
        }  // end of loop over integration point

        _bhe.template assembleRMatrices<ShapeFunction::NPOINTS>(
            idx_bhe_unknowns, matBHE_loc_R, _R_matrix, _R_pi_s_matrix,
            _R_s_matrix);
    }  // end of loop over BHE unknowns
}

template <typename ShapeFunction, typename BHEType>
void HeatTransportBHELocalAssemblerBHE<ShapeFunction, BHEType>::assemble(
    double const /*t*/, double const /*dt*/,
    std::vector<double> const& /*local_x*/,
    std::vector<double> const& /*local_x_prev*/,
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
    auto const& pipe_heat_conductions =
        _bhe.pipeHeatConductions(_section_index);
    auto const& pipe_advection_vectors =
        _bhe.pipeAdvectionVectors(_element_direction, _section_index);
    auto const& cross_section_areas = _bhe.crossSectionAreas(_section_index);

    // the mass and conductance matrix terms
    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        auto const& ip_data = _ip_data[ip];

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
}

template <typename ShapeFunction, typename BHEType>
void HeatTransportBHELocalAssemblerBHE<ShapeFunction, BHEType>::
    assembleWithJacobian(double const t, double const dt,
                         std::vector<double> const& local_x,
                         std::vector<double> const& local_x_prev,
                         std::vector<double>& local_rhs_data,
                         std::vector<double>& local_Jac_data)
{
    auto const local_matrix_size = local_x.size();
    // initialize x and x_prev
    auto x =
        Eigen::Map<BheLocalVectorType const>(local_x.data(), local_matrix_size);
    auto x_prev = Eigen::Map<BheLocalVectorType const>(local_x_prev.data(),
                                                       local_matrix_size);
    // initialize local_Jac and local_rhs
    auto local_Jac = MathLib::createZeroedMatrix<BheLocalMatrixType>(
        local_Jac_data, local_matrix_size, local_matrix_size);
    auto local_rhs = MathLib::createZeroedVector<BheLocalVectorType>(
        local_rhs_data, local_matrix_size);

    std::vector<double> local_M_data;
    std::vector<double> local_K_data;
    assemble(t, dt, local_x, local_x_prev, local_M_data, local_K_data,
             local_rhs_data /*not going to be used*/);

    // convert to matrix
    auto local_M = MathLib::toMatrix<BheLocalMatrixType>(
        local_M_data, local_matrix_size, local_matrix_size);
    auto local_K = MathLib::toMatrix<BheLocalMatrixType>(
        local_K_data, local_matrix_size, local_matrix_size);

    // Jac matrix and rhs vector operation
    local_Jac.noalias() += local_K + local_M / dt;
    local_rhs.noalias() -= local_K * x + local_M * (x - x_prev) / dt;

    local_M.setZero();
    local_K.setZero();
}

}  // namespace HeatTransportBHE
}  // namespace ProcessLib
