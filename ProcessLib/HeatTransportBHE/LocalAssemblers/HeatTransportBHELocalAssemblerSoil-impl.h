/**
 * \file
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <valarray>
#include <vector>

#include "MaterialLib/MPL/Medium.h"
#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "ProcessLib/HeatTransportBHE/HeatTransportBHEProcessData.h"
#include "ProcessLib/Utils/InitShapeMatrices.h"

#include "HeatTransportBHEProcessAssemblerInterface.h"
#include "SecondaryData.h"

namespace ProcessLib
{
namespace HeatTransportBHE
{
template <typename ShapeFunction, typename IntegrationMethod>
HeatTransportBHELocalAssemblerSoil<ShapeFunction, IntegrationMethod>::
    HeatTransportBHELocalAssemblerSoil(
        MeshLib::Element const& e,
        bool const is_axially_symmetric,
        unsigned const integration_order,
        HeatTransportBHEProcessData& process_data)
    : _process_data(process_data),
      _integration_method(integration_order),
      _element_id(e.getID())
{
    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    _ip_data.reserve(n_integration_points);
    _secondary_data.N.resize(n_integration_points);

    _shape_matrices = initShapeMatrices<ShapeFunction,
                                        ShapeMatricesType,
                                        IntegrationMethod,
                                        3 /* GlobalDim */>(
        e, is_axially_symmetric, _integration_method);

    ParameterLib::SpatialPosition x_position;
    x_position.setElementID(_element_id);

    // ip data initialization
    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        x_position.setIntegrationPoint(ip);

        // create the class IntegrationPointDataBHE in place
        auto const& sm = _shape_matrices[ip];
        double const w = _integration_method.getWeightedPoint(ip).getWeight() *
                         sm.integralMeasure * sm.detJ;
        _ip_data.push_back(
            {sm.N.transpose() * sm.N * w, sm.dNdx.transpose() * sm.dNdx * w});

        _secondary_data.N[ip] = sm.N;
    }
}

template <typename ShapeFunction, typename IntegrationMethod>
void HeatTransportBHELocalAssemblerSoil<ShapeFunction, IntegrationMethod>::
    assemble(double const t, double const /*dt*/,
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

    ParameterLib::SpatialPosition pos;
    pos.setElementID(_element_id);

    auto const& medium = *_process_data.media_map->getMedium(_element_id);
    auto const& solid_phase = medium.phase("Solid");
    auto const& liquid_phase = medium.phase("AqueousLiquid");

    MaterialPropertyLib::VariableArray vars;

    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        pos.setIntegrationPoint(ip);
        auto& ip_data = _ip_data[ip];
        auto const& M_ip = ip_data.NTN_product_times_w;
        auto const& K_ip = ip_data.dNdxTdNdx_product_times_w;

        auto const k_s =
            solid_phase
                .property(
                    MaterialPropertyLib::PropertyType::thermal_conductivity)
                .template value<double>(vars, pos, t);

        auto const heat_capacity_s =
            solid_phase
                .property(
                    MaterialPropertyLib::PropertyType::specific_heat_capacity)
                .template value<double>(vars, pos, t);

        auto const heat_capacity_f =
            liquid_phase
                .property(
                    MaterialPropertyLib::PropertyType::specific_heat_capacity)
                .template value<double>(vars, pos, t);

        auto const density_s =
            solid_phase.property(MaterialPropertyLib::PropertyType::density)
                .template value<double>(vars, pos, t);

        auto const density_f =
            liquid_phase.property(MaterialPropertyLib::PropertyType::density)
                .template value<double>(vars, pos, t);
        // for now only using the solid phase parameters

        auto const velocity =
            liquid_phase
                .property(MaterialPropertyLib::PropertyType::phase_velocity)
                .value(vars, pos, t);

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
