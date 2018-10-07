/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "HeatTransportBHELocalAssemblerSoil.h"

#include <valarray>
#include <vector>

#include <Eigen/Eigen>

#include "MathLib/LinAlg/Eigen/EigenMapTools.h"

#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "ProcessLib/Utils/InitShapeMatrices.h"

#include "ProcessLib/HeatTransportBHE/HeatTransportBHEProcessData.h"

// #include "IntegrationPointDataMatrix.h"
#include "HeatTransportBHEProcessAssemblerInterface.h"
#include "SecondaryData.h"

namespace ProcessLib
{
namespace HeatTransportBHE
{
template <typename ShapeFunction, typename IntegrationMethod, int GlobalDim>
HeatTransportBHELocalAssemblerSoil<ShapeFunction,
                                   IntegrationMethod,
                                   GlobalDim>::
    HeatTransportBHELocalAssemblerSoil(
        MeshLib::Element const& e,
        std::vector<unsigned> const& dofIndex_to_localIndex,
        bool const is_axially_symmetric,
        unsigned const integration_order,
        HeatTransportBHEProcessData& process_data)
    : HeatTransportBHELocalAssemblerInterface(
          ShapeFunction::NPOINTS * GlobalDim, dofIndex_to_localIndex),
      _process_data(process_data),
      _integration_method(integration_order),
      element_id(e.getID()),
      _is_axially_symmetric(is_axially_symmetric)
{
    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    _ip_data.reserve(n_integration_points);
    _secondary_data.N.resize(n_integration_points);

    _shape_matrices = initShapeMatrices<ShapeFunction,
                                        ShapeMatricesType,
                                        IntegrationMethod,
                                        GlobalDim>(
        e, is_axially_symmetric, _integration_method);

    SpatialPosition x_position;
    x_position.setElementID(element_id);

    // ip data initialization
    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        x_position.setIntegrationPoint(ip);

        // create the class IntegrationPointDataBHE in place
        _ip_data.emplace_back();
        auto const& sm = _shape_matrices[ip];
        auto& ip_data = _ip_data[ip];
        double const w = _integration_method.getWeightedPoint(ip).getWeight() *
                         sm.integralMeasure * sm.detJ;
        ip_data.NTN_product_times_w = sm.N.transpose() * sm.N * w;
        ip_data.dNdxTdNdx_product_times_w = sm.dNdx.transpose() * sm.dNdx * w;

        _secondary_data.N[ip] = sm.N;
    }
}

template <typename ShapeFunction, typename IntegrationMethod, int GlobalDim>
void HeatTransportBHELocalAssemblerSoil<
    ShapeFunction,
    IntegrationMethod,
    GlobalDim>::assemble(double const t,
                         std::vector<double> const& local_x,
                         std::vector<double>& local_M_data,
                         std::vector<double>& local_K_data,
                         std::vector<double>& /*local_b_data*/)
{
    assert(local_x.size() == ShapeFunction::NPOINTS);
    (void)local_x;  // Avoid unused arg warning.

    auto local_M = MathLib::createZeroedMatrix<NodalMatrixType>(
        local_M_data, ShapeFunction::NPOINTS, ShapeFunction::NPOINTS);
    auto local_K = MathLib::createZeroedMatrix<NodalMatrixType>(
        local_K_data, ShapeFunction::NPOINTS, ShapeFunction::NPOINTS);

    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    SpatialPosition pos;
    pos.setElementID(element_id);

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        pos.setIntegrationPoint(ip);
        auto& ip_data = _ip_data[ip];
        auto const& M_ip = ip_data.NTN_product_times_w;
        auto const& K_ip = ip_data.dNdxTdNdx_product_times_w;

        // auto const k_f = _process_data.thermal_conductivity_fluid(t, pos)[0];
        // auto const k_g = _process_data.thermal_conductivity_gas(t, pos)[0];
        auto const k_s = _process_data.thermal_conductivity_solid(t, pos)[0];

        // auto const heat_capacity_f = _process_data.heat_capacity_fluid(t,
        // pos)[0]; auto const heat_capacity_g =
        // _process_data.heat_capacity_gas(t, pos)[0];
        auto const heat_capacity_s =
            _process_data.heat_capacity_solid(t, pos)[0];

        // auto const density_f = _process_data.density_fluid(t, pos)[0];
        // auto const density_g = _process_data.density_gas(t, pos)[0];
        auto const density_s = _process_data.density_solid(t, pos)[0];

        // for now only using the solid phase parameters

        // assemble Conductance matrix
        local_K.noalias() += K_ip * k_s;

        // assemble Mass matrix
        local_M.noalias() += M_ip * density_s * heat_capacity_s;
    }

    // debugging
    // std::string sep = "\n----------------------------------------\n";
    // Eigen::IOFormat CleanFmt(4, 0, ", ", "\n", "[", "]");
    // std::cout << local_K.format(CleanFmt) << sep;
    // std::cout << local_M.format(CleanFmt) << sep;
}
}  // namespace HeatTransportBHE
}  // namespace ProcessLib
