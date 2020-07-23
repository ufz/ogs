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

#include <valarray>
#include <vector>

#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/MPL/Utils/FormEigenTensor.h"
#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "NumLib/Function/Interpolation.h"
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

    _shape_matrices =
        initShapeMatrices<ShapeFunction, ShapeMatricesType, 3 /* GlobalDim */>(
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
        _ip_data.push_back({sm.N, sm.dNdx, w});

        _secondary_data.N[ip] = sm.N;
    }
}

template <typename ShapeFunction, typename IntegrationMethod>
void HeatTransportBHELocalAssemblerSoil<ShapeFunction, IntegrationMethod>::
    assemble(double const t, double const dt,
             std::vector<double> const& local_x,
             std::vector<double> const& /*local_xdot*/,
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
        auto const& N = ip_data.N;
        auto const& dNdx = ip_data.dNdx;
        auto const& w = ip_data.integration_weight;

        double T_int_pt = 0.0;
        NumLib::shapeFunctionInterpolate(local_x, N, T_int_pt);

        vars[static_cast<int>(MaterialPropertyLib::Variable::temperature)] =
            T_int_pt;

        // for now only using the solid and liquid phase parameters
        auto const density_s =
            solid_phase.property(MaterialPropertyLib::PropertyType::density)
                .template value<double>(vars, pos, t, dt);

        auto const heat_capacity_s =
            solid_phase
                .property(
                    MaterialPropertyLib::PropertyType::specific_heat_capacity)
                .template value<double>(vars, pos, t, dt);

        auto const density_f =
            liquid_phase.property(MaterialPropertyLib::PropertyType::density)
                .template value<double>(vars, pos, t, dt);

        auto const heat_capacity_f =
            liquid_phase
                .property(
                    MaterialPropertyLib::PropertyType::specific_heat_capacity)
                .template value<double>(vars, pos, t, dt);

        auto const porosity =
            medium.property(MaterialPropertyLib::PropertyType::porosity)
                .template value<double>(vars, pos, t, dt);

        auto const velocity =
            liquid_phase
                .property(MaterialPropertyLib::PropertyType::phase_velocity)
                .template value<Eigen::Vector3d>(vars, pos, t, dt);

        // calculate the hydrodynamic thermodispersion tensor
        auto const thermal_conductivity =
            MaterialPropertyLib::formEigenTensor<3>(
                medium
                    .property(
                        MaterialPropertyLib::PropertyType::thermal_conductivity)
                    .value(vars, pos, t, dt));

        auto thermal_conductivity_dispersivity = thermal_conductivity;

        double const velocity_magnitude = velocity.norm();

        if (velocity_magnitude >= std::numeric_limits<double>::epsilon())
        {
            auto const thermal_dispersivity_longitudinal =
                medium
                    .property(MaterialPropertyLib::PropertyType::
                                  thermal_longitudinal_dispersivity)
                    .template value<double>();
            auto const thermal_dispersivity_transversal =
                medium
                    .property(MaterialPropertyLib::PropertyType::
                                  thermal_transversal_dispersivity)
                    .template value<double>();

            auto const thermal_dispersivity =
                density_f * heat_capacity_f *
                (thermal_dispersivity_transversal * velocity_magnitude *
                     Eigen::Matrix3d::Identity() +
                 (thermal_dispersivity_longitudinal -
                  thermal_dispersivity_transversal) /
                     velocity_magnitude * velocity * velocity.transpose());
            thermal_conductivity_dispersivity += thermal_dispersivity;
        }

        // assemble Conductance matrix
        local_K.noalias() +=
            (dNdx.transpose() * thermal_conductivity_dispersivity * dNdx +
             N.transpose() * velocity.transpose() * dNdx * density_f *
                 heat_capacity_f) *
            w;

        // assemble Mass matrix
        local_M.noalias() += N.transpose() * N * w *
                             (density_s * heat_capacity_s * (1 - porosity) +
                              density_f * heat_capacity_f * porosity);
    }

    // debugging
    // std::string sep = "\n----------------------------------------\n";
    // Eigen::IOFormat CleanFmt(4, 0, ", ", "\n", "[", "]");
    // std::cout << local_K.format(CleanFmt) << sep;
    // std::cout << local_M.format(CleanFmt) << sep;
}
}  // namespace HeatTransportBHE
}  // namespace ProcessLib
