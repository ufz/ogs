/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <vector>

#include "HeatTransportBHEProcessAssemblerInterface.h"
#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/MPL/Utils/FormEigenTensor.h"
#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "NumLib/Fem/InitShapeMatrices.h"
#include "NumLib/Fem/Interpolation.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "ProcessLib/HeatTransportBHE/HeatTransportBHEProcess.h"
#include "ProcessLib/HeatTransportBHE/HeatTransportBHEProcessData.h"
#include "SecondaryData.h"

namespace ProcessLib
{
namespace HeatTransportBHE
{
template <typename ShapeFunction>
HeatTransportBHELocalAssemblerSoil<ShapeFunction>::
    HeatTransportBHELocalAssemblerSoil(
        MeshLib::Element const& e,
        NumLib::GenericIntegrationMethod const& integration_method,
        bool const is_axially_symmetric,
        HeatTransportBHEProcessData& process_data)
    : _process_data(process_data),
      _integration_method(integration_method),
      _element(e)
{
    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    _ip_data.reserve(n_integration_points);
    _secondary_data.N.resize(n_integration_points);

    _shape_matrices =
        NumLib::initShapeMatrices<ShapeFunction, ShapeMatricesType,
                                  3 /* GlobalDim */>(e, is_axially_symmetric,
                                                     _integration_method);

    // ip data initialization
    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        // create the class IntegrationPointDataBHE in place
        auto const& sm = _shape_matrices[ip];
        double const w = _integration_method.getWeightedPoint(ip).getWeight() *
                         sm.integralMeasure * sm.detJ;
        _ip_data.push_back({sm.N, sm.dNdx, w});

        _secondary_data.N[ip] = sm.N;
    }
}

template <typename ShapeFunction>
void HeatTransportBHELocalAssemblerSoil<ShapeFunction>::assemble(
    double const t, double const dt, std::vector<double> const& local_x,
    std::vector<double> const& /*local_x_prev*/,
    std::vector<double>& local_M_data, std::vector<double>& local_K_data,
    std::vector<double>& /*local_b_data*/)
{
    assert(local_x.size() == ShapeFunction::NPOINTS);
    (void)local_x;  // Avoid unused arg warning.

    auto local_M = MathLib::createZeroedMatrix<NodalMatrixType>(
        local_M_data, ShapeFunction::NPOINTS, ShapeFunction::NPOINTS);
    auto local_K = MathLib::createZeroedMatrix<NodalMatrixType>(
        local_K_data, ShapeFunction::NPOINTS, ShapeFunction::NPOINTS);

    auto const& medium = *_process_data.media_map.getMedium(_element.getID());
    auto const& solid_phase = medium.phase("Solid");
    auto const& liquid_phase = medium.phase("AqueousLiquid");

    MaterialPropertyLib::VariableArray vars;

    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        auto& ip_data = _ip_data[ip];
        auto const& N = ip_data.N;
        auto const& dNdx = ip_data.dNdx;
        auto const& w = ip_data.integration_weight;

        ParameterLib::SpatialPosition const pos{
            std::nullopt, _element.getID(),
            MathLib::Point3d(NumLib::interpolateCoordinates<ShapeFunction,
                                                            ShapeMatricesType>(
                _element, N))};

        double T_int_pt = 0.0;
        NumLib::shapeFunctionInterpolate(local_x, N, T_int_pt);

        vars.temperature = T_int_pt;

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

    if (_process_data._mass_lumping)
    {
        // only mass lumping at the BHE connected soil elements
        if (_process_data.mass_lumping_soil_elements[_element.getID()])
        {
            local_M = local_M.colwise().sum().eval().asDiagonal();
        }
    }

    //  debugging
    //  std::string sep = "\n----------------------------------------\n";
    //  Eigen::IOFormat CleanFmt(4, 0, ", ", "\n", "[", "]");
    //  std::cout << local_K.format(CleanFmt) << sep;
    //  std::cout << local_M.format(CleanFmt) << sep;
}

template <typename ShapeFunction>
void HeatTransportBHELocalAssemblerSoil<ShapeFunction>::assembleWithJacobian(
    double const t, double const dt, std::vector<double> const& local_x,
    std::vector<double> const& local_x_prev,
    std::vector<double>& local_rhs_data, std::vector<double>& local_Jac_data)
{
    assert(local_x.size() == ShapeFunction::NPOINTS);
    auto const local_matrix_size = local_x.size();
    // initialize x and x_prev
    auto x =
        Eigen::Map<NodalVectorType const>(local_x.data(), local_matrix_size);
    auto x_prev = Eigen::Map<NodalVectorType const>(local_x_prev.data(),
                                                    local_matrix_size);
    // initialize local_Jac and local_rhs
    auto local_Jac = MathLib::createZeroedMatrix<NodalMatrixType>(
        local_Jac_data, local_matrix_size, local_matrix_size);
    auto local_rhs = MathLib::createZeroedVector<NodalVectorType>(
        local_rhs_data, local_matrix_size);

    std::vector<double> local_M_data;
    std::vector<double> local_K_data;
    assemble(t, dt, local_x, local_x_prev, local_M_data, local_K_data,
             local_rhs_data /*not going to be used*/);

    // convert to matrix
    auto local_M = MathLib::toMatrix<NodalMatrixType>(
        local_M_data, local_matrix_size, local_matrix_size);
    auto local_K = MathLib::toMatrix<NodalMatrixType>(
        local_K_data, local_matrix_size, local_matrix_size);

    // Jac matrix and rhs vector operation
    local_Jac.noalias() += local_K + local_M / dt;
    local_rhs.noalias() -= local_K * x + local_M * (x - x_prev) / dt;

    local_M.setZero();
    local_K.setZero();
}

}  // namespace HeatTransportBHE
}  // namespace ProcessLib
