// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <spdlog/fmt/fmt.h>

#include <algorithm>
#include <cmath>
#include <limits>

#include "GenericNaturalBoundaryConditionLocalAssembler.h"
#include "MaterialLib/MPL/MaterialSpatialDistributionMap.h"
#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/MPL/Utils/DriftFluxModel.h"
#include "MeshLib/PropertyVector.h"
#include "NumLib/DOF/DOFTableUtil.h"
#include "NumLib/Exceptions.h"
#include "NumLib/Fem/Interpolation.h"
#include "NumLib/IndexValueVector.h"
#include "ParameterLib/MeshNodeParameter.h"

namespace ProcessLib
{

struct WellboreCompensateCoefficients
{
    double pressure;
    double velocity;
    double enthalpy;
};

struct WellboreCompensateNeumannBoundaryConditionData
{
    WellboreCompensateCoefficients coefficients;

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
        auto const& liquid_phase =
            medium.phase(MaterialPropertyLib::PhaseName::AqueousLiquid);
        auto const& gas_phase =
            medium.phase(MaterialPropertyLib::PhaseName::Gas);

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

            double const dryness = std::clamp(
                (enthalpy_int_pt - h_sat_liq_w) / (h_sat_vap_w - h_sat_liq_w),
                0., 1.);

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
            double const C_0 =
                MaterialPropertyLib::driftFluxProfileParameter(dryness);

            // drift flux velocity
            double const u_gu = MaterialPropertyLib::driftFluxVelocity(
                dryness, T_int_pt, vapour_water_density, liquid_water_density);

            MaterialPropertyLib::DriftFluxState const drift_flux_state{
                .dryness = dryness,
                .vapour_water_density = vapour_water_density,
                .liquid_water_density = liquid_water_density,
                .v_mix = velocity_int_pt,
                .C_0 = C_0,
                .u_gu = u_gu};

            // solving void fraction of vapour: Rouhani-Axelsson
            auto const alpha_solution =
                MaterialPropertyLib::computeVapourVoidFraction(
                    drift_flux_state);

            if (!alpha_solution)
            {
                throw NumLib::AssemblyException(fmt::format(
                    "The drift-flux closure of the WellboreCompensateNeumann "
                    "boundary condition has no admissible vapour void fraction "
                    "in element {:d}, integration point {:d}: pressure {:g}, "
                    "mixture velocity {:g}, specific enthalpy {:g}, "
                    "temperature {:g}, {}",
                    _element.getID(), ip, pressure_int_pt, velocity_int_pt,
                    enthalpy_int_pt, T_int_pt,
                    MaterialPropertyLib::voidFractionClosureDiagnostics(
                        drift_flux_state)));
            }

            double const alpha = *alpha_solution;

            if (alpha == 0)
            {
                liquid_water_density =
                    liquid_phase
                        .property(MaterialPropertyLib::PropertyType::density)
                        .template value<double>(vars, pos, 0, 0);
            }

            double const mix_density = vapour_water_density * alpha +
                                       liquid_water_density * (1 - alpha);

            double const gamma = MaterialPropertyLib::mixtureSlipParameter(
                alpha, drift_flux_state);

            double const neumann_ip_values =
                _data.coefficients.pressure * mix_density * velocity_int_pt +
                _data.coefficients.velocity *
                    (mix_density * velocity_int_pt * velocity_int_pt + gamma) +
                _data.coefficients.enthalpy * mix_density * velocity_int_pt *
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
