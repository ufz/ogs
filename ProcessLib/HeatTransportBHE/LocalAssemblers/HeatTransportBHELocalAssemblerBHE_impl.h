/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <Eigen/Eigen>
#include "HeatTransportBHELocalAssemblerBHE.h"
#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "ProcessLib/HeatTransportBHE/LocalAssemblers/IntegrationPointDataBHE.h"
#include "ProcessLib/Utils/InitShapeMatrices.h"

namespace ProcessLib
{
namespace HeatTransportBHE
{
using namespace BHE;

template <typename ShapeFunction, typename IntegrationMethod, int BHE_Dim>
HeatTransportBHELocalAssemblerBHE<ShapeFunction, IntegrationMethod, BHE_Dim>::
    HeatTransportBHELocalAssemblerBHE(
        MeshLib::Element const& e,
        std::size_t const /*local_matrix_size*/,
        std::vector<unsigned> const& dofIndex_to_localIndex,
        bool const is_axially_symmetric,
        unsigned const integration_order,
        HeatTransportBHEProcessData& process_data)
    : HeatTransportBHELocalAssemblerInterface(
          ShapeFunction::NPOINTS * BHE_Dim,  // no intersection
          dofIndex_to_localIndex),
      _process_data(process_data),
      _integration_method(integration_order),
      element_id(e.getID())
{
    // need to make sure that the BHE elements are one-dimensional
    assert(e.getDimension() == 1);

    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    _ip_data.reserve(n_integration_points);
    _secondary_data.N.resize(n_integration_points);

    auto const shape_matrices =
        initShapeMatrices<ShapeFunction, ShapeMatricesType, IntegrationMethod,
                          BHE_Dim>(e, is_axially_symmetric,
                                   _integration_method);

    auto mat_id = (*_process_data._mesh_prop_materialIDs)[e.getID()];
    auto BHE_id = _process_data._map_materialID_to_BHE_ID[mat_id];

    SpatialPosition x_position;
    x_position.setElementID(element_id);

    // ip data initialization
    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        x_position.setIntegrationPoint(ip);

        IntegrationPointDataBHE<ShapeMatricesType> int_Point_Data_BHE(
            *(_process_data._vec_BHE_property[BHE_id]));
        _ip_data.emplace_back(int_Point_Data_BHE);
        auto const& sm = shape_matrices[ip];
        auto& ip_data = _ip_data[ip];
        ip_data.integration_weight =
            _integration_method.getWeightedPoint(ip).getWeight() *
            sm.integralMeasure * sm.detJ;
        ip_data.N = sm.N;
        ip_data.dNdx = sm.dNdx;

        _secondary_data.N[ip] = sm.N;
    }

    const int BHE_n_unknowns = _ip_data[0]._bhe_instance.getNumUnknowns();
    const int NumRMatrixRows = ShapeFunction::NPOINTS * BHE_n_unknowns;
    _R_matrix.setZero(NumRMatrixRows, NumRMatrixRows);
    _R_pi_s_matrix.setZero(NumRMatrixRows, ShapeFunction::NPOINTS);
    _R_s_matrix.setZero(ShapeFunction::NPOINTS, ShapeFunction::NPOINTS);
    // formulate the local BHE R matrix
    Eigen::MatrixXd matBHE_loc_R =
        Eigen::MatrixXd::Zero(ShapeFunction::NPOINTS, ShapeFunction::NPOINTS);
    for (int idx_bhe_unknowns = 0; idx_bhe_unknowns < BHE_n_unknowns;
         idx_bhe_unknowns++)
    {
        matBHE_loc_R.setZero();
        // Loop over Gauss points
        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            x_position.setIntegrationPoint(ip);

            auto const& N = _ip_data[ip].N;
            auto const& w = _ip_data[ip].integration_weight;

            // get coefficient of R matrix for corresponding BHE.
            auto R_coeff = _process_data._vec_BHE_property[BHE_id]
                               ->getBoundaryHeatExchangeCoeff(idx_bhe_unknowns);

            // calculate mass matrix for current unknown
            matBHE_loc_R += N.transpose() * R_coeff * N * w;
        }  // end of loop over integration point

        // The following assembly action is according to Diersch (2013) FEFLOW
        // book please refer to M.127 and M.128 on page 955 and 956
        _process_data._vec_BHE_property[BHE_id]->setRMatrices(
            idx_bhe_unknowns, ShapeFunction::NPOINTS, matBHE_loc_R, _R_matrix,
            _R_pi_s_matrix, _R_s_matrix);

    }  // end of loop over BHE unknowns

    // debugging
    // std::string sep = "\n----------------------------------------\n";
    // Eigen::IOFormat CleanFmt(4, 0, ", ", "\n", "[", "]");
    // std::cout << "_R_matrix: \n" << sep;
    // std::cout << _R_matrix.format(CleanFmt) << sep;
    // std::cout << "_R_s_matrix: \n" << sep;
    // std::cout << _R_s_matrix.format(CleanFmt) << sep;
    // std::cout << "_R_pi_s_matrix: \n" << sep;
    // std::cout << _R_pi_s_matrix.format(CleanFmt) << sep;
}

template <typename ShapeFunction, typename IntegrationMethod, int BHE_Dim>
void HeatTransportBHELocalAssemblerBHE<ShapeFunction, IntegrationMethod,
                                       BHE_Dim>::
    assemble(
        double const /*t*/, std::vector<double> const& local_x,
        std::vector<double>& local_M_data, std::vector<double>& local_K_data,
        std::vector<double>& /*local_b_data*/)  // local b vector is not touched
{
    auto const local_matrix_size = local_x.size();
    const int BHE_n_unknowns = _ip_data[0]._bhe_instance.getNumUnknowns();
    // plus one because the soil temperature is included in local_x
    assert(local_matrix_size == ShapeFunction::NPOINTS * (BHE_n_unknowns + 1));

    auto local_M = MathLib::createZeroedMatrix<BheLocalMatrixType>(
        local_M_data, local_matrix_size, local_matrix_size);
    auto local_K = MathLib::createZeroedMatrix<BheLocalMatrixType>(
        local_K_data, local_matrix_size, local_matrix_size);

    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    SpatialPosition x_position;
    x_position.setElementID(element_id);

    int shift = 0;
    static const int local_BHE_matrix_size =
        ShapeFunction::NPOINTS * BHE_n_unknowns;
    static const int shift_start = ShapeFunction::NPOINTS;

    // the mass and conductance matrix terms
    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        x_position.setIntegrationPoint(ip);
        auto& ip_data = _ip_data[ip];

        auto const& w = ip_data.integration_weight;
        auto const& N = ip_data.N;
        auto const& dNdx = ip_data.dNdx;

        // looping over all unknowns.
        for (int idx_bhe_unknowns = 0; idx_bhe_unknowns < BHE_n_unknowns;
             idx_bhe_unknowns++)
        {
            // get coefficient of mass from corresponding BHE.
            auto& mass_coeff = ip_data._vec_mass_coefficients[idx_bhe_unknowns];
            auto& laplace_mat = ip_data._vec_mat_Laplace[idx_bhe_unknowns];
            auto& advection_vec =
                ip_data._vec_Advection_vectors[idx_bhe_unknowns];

            // calculate shift.
            shift = shift_start + ShapeFunction::NPOINTS * idx_bhe_unknowns;
            // local M
            local_M
                .template block<ShapeFunction::NPOINTS, ShapeFunction::NPOINTS>(
                    shift, shift)
                .noalias() += N.transpose() * mass_coeff * N * w;

            // local K
            // laplace part
            local_K
                .block(shift, shift, ShapeFunction::NPOINTS,
                       ShapeFunction::NPOINTS)
                .noalias() += dNdx.transpose() * laplace_mat * dNdx * w;
            // advection part
            local_K
                .block(shift, shift, ShapeFunction::NPOINTS,
                       ShapeFunction::NPOINTS)
                .noalias() +=
                N.transpose() * advection_vec.transpose() * dNdx * w;
        }
    }

    // add the R matrix to local_K
    local_K.block(shift_start, shift_start, local_BHE_matrix_size,
                  local_BHE_matrix_size) += _R_matrix;

    // add the R_pi_s matrix to local_K
    local_K.block(shift_start, 0, local_BHE_matrix_size, shift_start) +=
        _R_pi_s_matrix;
    local_K.block(0, shift_start, shift_start, local_BHE_matrix_size) +=
        _R_pi_s_matrix.transpose();

    // add the R_s matrix to local_K
    local_K.block(0, 0, shift_start, shift_start) += 2.0 * _R_s_matrix;

    // debugging
    // std::string sep =
    //     "\n----------------------------------------\n";
    // Eigen::IOFormat CleanFmt(4, 0, ", ", "\n", "[", "]");
    // std::cout << local_K.format(CleanFmt) << sep;
    // std::cout << local_M.format(CleanFmt) << sep;
}

template <typename ShapeFunction, typename IntegrationMethod, int BHE_Dim>
void HeatTransportBHELocalAssemblerBHE<
    ShapeFunction, IntegrationMethod,
    BHE_Dim>::postTimestepConcrete(std::vector<double> const& /*local_x*/)
{
}
}  // namespace HeatTransportBHE
}  // namespace ProcessLib
