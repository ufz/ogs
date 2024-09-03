/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <Eigen/LU>
#include <limits>

#include "GenericNaturalBoundaryConditionLocalAssembler.h"
#include "MaterialLib/MPL/MaterialSpatialDistributionMap.h"
#include "MaterialLib/MPL/Medium.h"
#include "MeshLib/PropertyVector.h"
#include "NumLib/DOF/DOFTableUtil.h"
#include "NumLib/Fem/Interpolation.h"
#include "NumLib/IndexValueVector.h"
#include "NumLib/NewtonRaphson.h"
#include "ParameterLib/MeshNodeParameter.h"

namespace ProcessLib
{
struct WellboreCompensateNeumannBoundaryConditionData
{
    double const pressure_coefficient;
    double const velocity_coefficient;
    double const enthalpy_coefficient;

    // Used for mapping boundary nodes to bulk nodes.
    std::unique_ptr<NumLib::LocalToGlobalIndexMap> dof_table_boundary_pressure;
    std::unique_ptr<NumLib::LocalToGlobalIndexMap> dof_table_boundary_velocity;
    std::unique_ptr<NumLib::LocalToGlobalIndexMap> dof_table_boundary_enthalpy;

    MaterialPropertyLib::MaterialSpatialDistributionMap media_map;
};

template <typename ShapeFunction, int GlobalDim>
class WellboreCompensateNeumannBoundaryConditionLocalAssembler final
    : public GenericNaturalBoundaryConditionLocalAssembler<ShapeFunction,
                                                           GlobalDim>
{
    using Base =
        GenericNaturalBoundaryConditionLocalAssembler<ShapeFunction, GlobalDim>;
    using NodalVectorType = typename Base::NodalVectorType;
    using NodalMatrixType = typename Base::NodalMatrixType;

public:
    /// The neumann_bc_term factor is directly integrated into the local
    /// element matrix.
    WellboreCompensateNeumannBoundaryConditionLocalAssembler(
        MeshLib::Element const& e,
        std::size_t const local_matrix_size,
        NumLib::GenericIntegrationMethod const& integration_method,
        bool const is_axially_symmetric,
        WellboreCompensateNeumannBoundaryConditionData const& data)
        : Base(e, is_axially_symmetric, integration_method),
          _element(e),
          _data(data),
          _local_matrix_size(local_matrix_size)
    {
    }

    void assemble(std::size_t const mesh_item_id,
                  NumLib::LocalToGlobalIndexMap const& dof_table_boundary,
                  double const /*t*/, std::vector<GlobalVector*> const& x,
                  int const process_id, GlobalMatrix* /*K*/, GlobalVector& b,
                  GlobalMatrix* /*Jac*/) override
    {
        NodalVectorType _local_rhs(_local_matrix_size);
        _local_rhs.setZero();

        unsigned const n_integration_points =
            Base::_integration_method.getNumberOfPoints();

        auto const indices_current_variable =
            NumLib::getIndices(mesh_item_id, dof_table_boundary);
        auto const indices_pressure = NumLib::getIndices(
            mesh_item_id, *_data.dof_table_boundary_pressure);
        auto const indices_velocity = NumLib::getIndices(
            mesh_item_id, *_data.dof_table_boundary_velocity);
        auto const indices_enthalpy = NumLib::getIndices(
            mesh_item_id, *_data.dof_table_boundary_enthalpy);

        std::vector<double> const local_pressure =
            x[process_id]->get(indices_pressure);
        std::vector<double> const local_velocity =
            x[process_id]->get(indices_velocity);
        std::vector<double> const local_enthalpy =
            x[process_id]->get(indices_enthalpy);

        auto const& medium = *_data.media_map.getMedium(_element.getID());
        auto const& liquid_phase = medium.phase("AqueousLiquid");
        auto const& gas_phase = medium.phase("Gas");

        ParameterLib::SpatialPosition pos;
        pos.setElementID(_element.getID());

        MaterialPropertyLib::VariableArray vars;

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            auto const& n_and_weight = Base::_ns_and_weights[ip];
            auto const& N = n_and_weight.N;
            auto const& w = n_and_weight.weight;

            double pressure_int_pt = 0.0;
            double velocity_int_pt = 0.0;
            double enthalpy_int_pt = 0.0;

            NumLib::shapeFunctionInterpolate(local_pressure, N,
                                             pressure_int_pt);
            NumLib::shapeFunctionInterpolate(local_velocity, N,
                                             velocity_int_pt);
            NumLib::shapeFunctionInterpolate(local_enthalpy, N,
                                             enthalpy_int_pt);

            vars.liquid_phase_pressure = pressure_int_pt;
            vars.enthalpy = enthalpy_int_pt;

            double liquid_water_density =
                liquid_phase
                    .property(
                        MaterialPropertyLib::PropertyType::saturation_density)
                    .template value<double>(vars, pos, 0, 0);

            double const vapour_water_density =
                gas_phase
                    .property(
                        MaterialPropertyLib::PropertyType::saturation_density)
                    .template value<double>(vars, pos, 0, 0);

            double const h_sat_liq_w =
                liquid_phase
                    .property(
                        MaterialPropertyLib::PropertyType::saturation_enthalpy)
                    .template value<double>(vars, pos, 0, 0);

            double const h_sat_vap_w =
                gas_phase
                    .property(
                        MaterialPropertyLib::PropertyType::saturation_enthalpy)
                    .template value<double>(vars, pos, 0, 0);

            double const dryness =
                std::max(0., (enthalpy_int_pt - h_sat_liq_w) /
                                 (h_sat_vap_w - h_sat_liq_w));

            double const T_int_pt =
                (dryness == 0)
                    ? liquid_phase
                          .property(
                              MaterialPropertyLib::PropertyType::temperature)
                          .template value<double>(vars, pos, 0, 0)
                    : gas_phase
                          .property(MaterialPropertyLib::PropertyType::
                                        saturation_temperature)
                          .template value<double>(vars, pos, 0, 0);

            vars.temperature = T_int_pt;

            // For the calculation of the void fraction of vapour,
            // see Rohuani, Z., and E. Axelsson. "Calculation of volume void
            // fraction in a subcooled and quality region." International
            // Journal of Heat and Mass Transfer 17 (1970): 383-393.

            // Profile parameter of drift flux
            double const C_0 = 1 + 0.12 * (1 - dryness);

            // For the surface tension calculation, see
            // Cooper, J. R., and R. B. Dooley. "IAPWS release on surface
            // tension of ordinary water substance." International Association
            // for the Properties of Water and Steam (1994).
            double const sigma_gl = 0.2358 *
                                    std::pow((1 - T_int_pt / 647.096), 1.256) *
                                    (1 - 0.625 * (1 - T_int_pt / 647.096));
            // drift flux velocity
            double const u_gu =
                1.18 * (1 - dryness) *
                std::pow((9.81) * sigma_gl *
                             (liquid_water_density - vapour_water_density),
                         0.25) /
                std::pow(liquid_water_density, 0.5);

            // solving void fraction of vapour using local Newton
            // iteration.
            double alpha = 0;
            if (dryness != 0)
            {
                // Local Newton solver
                using LocalJacobianMatrix =
                    Eigen::Matrix<double, 1, 1, Eigen::RowMajor>;
                using LocalResidualVector = Eigen::Matrix<double, 1, 1>;
                using LocalUnknownVector = Eigen::Matrix<double, 1, 1>;
                LocalJacobianMatrix J_loc;

                Eigen::PartialPivLU<LocalJacobianMatrix> linear_solver(1);

                auto const update_residual = [&](LocalResidualVector& residual)
                {
                    residual(0) = dryness * liquid_water_density *
                                      (alpha * vapour_water_density +
                                       (1 - alpha) * liquid_water_density) *
                                      velocity_int_pt -
                                  alpha * C_0 * dryness * liquid_water_density *
                                      (alpha * vapour_water_density +
                                       (1 - alpha) * liquid_water_density) *
                                      velocity_int_pt -
                                  alpha * C_0 * (1 - dryness) *
                                      vapour_water_density *
                                      (alpha * vapour_water_density +
                                       (1 - alpha) * liquid_water_density) *
                                      velocity_int_pt -
                                  alpha * vapour_water_density *
                                      liquid_water_density * u_gu;
                };

                auto const update_jacobian = [&](LocalJacobianMatrix& jacobian)
                {
                    jacobian(0) =
                        dryness * liquid_water_density * velocity_int_pt *
                            (vapour_water_density - liquid_water_density) -
                        (C_0 * dryness * liquid_water_density +
                         C_0 * (1 - dryness) * vapour_water_density) *
                            (2 * alpha * vapour_water_density +
                             (1 - 2 * alpha) * liquid_water_density) *
                            velocity_int_pt -
                        vapour_water_density * liquid_water_density * u_gu;
                };

                auto const update_solution =
                    [&](LocalUnknownVector const& increment)
                {
                    // increment solution vectors
                    alpha += increment[0];
                };

                const int maximum_iterations(20);
                const double residuum_tolerance(1.e-10);
                const double increment_tolerance(0);

                auto newton_solver = NumLib::NewtonRaphson(
                    linear_solver, update_jacobian, update_residual,
                    update_solution,
                    {maximum_iterations, residuum_tolerance,
                     increment_tolerance});

                auto const success_iterations = newton_solver.solve(J_loc);

                if (!success_iterations)
                {
                    WARN(
                        "Attention! Steam void fraction has not been correctly "
                        "calculated!");
                }
            }

            if (alpha == 0)
            {
                liquid_water_density =
                    liquid_phase
                        .property(MaterialPropertyLib::PropertyType::density)
                        .template value<double>(vars, pos, 0, 0);
            }

            double const mix_density = vapour_water_density * alpha +
                                       liquid_water_density * (1 - alpha);

            // slip parameter between two phases
            double const gamma =
                alpha * liquid_water_density * vapour_water_density *
                mix_density / (1 - alpha) /
                std::pow((alpha * C_0 * vapour_water_density +
                          (1 - alpha * C_0) * liquid_water_density),
                         2) *
                std::pow((C_0 - 1) * velocity_int_pt + u_gu, 2);

            double const neumann_ip_values =
                _data.pressure_coefficient * mix_density * velocity_int_pt +
                _data.velocity_coefficient *
                    (mix_density * velocity_int_pt * velocity_int_pt + gamma) +
                _data.enthalpy_coefficient * mix_density * velocity_int_pt *
                    velocity_int_pt * velocity_int_pt * 0.5;
            _local_rhs.noalias() += N.transpose() * neumann_ip_values * w;
        }

        b.add(indices_current_variable, _local_rhs);
    }

private:
    MeshLib::Element const& _element;
    WellboreCompensateNeumannBoundaryConditionData const& _data;
    typename Base::NodalVectorType _local_matrix_size;
};

}  // namespace ProcessLib
