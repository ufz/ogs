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

#include "BaseLib/Error.h"

#include "BHECommon.h"
#include "FlowAndTemperatureControl.h"
#include "PipeConfiguration1U.h"

namespace ProcessLib
{
namespace HeatTransportBHE
{
namespace BHE
{
class BHE_1U final : public BHECommon
{
public:
    BHE_1U(BoreholeGeometry const& borehole,
           RefrigerantProperties const& refrigerant,
           GroutParameters const& grout,
           FlowAndTemperatureControl const& flowAndTemperatureControl,
           PipeConfiguration1U const& pipes)
        : BHECommon{borehole, refrigerant, grout, flowAndTemperatureControl},
          _pipes(pipes)
    {
        // Initialize thermal resistances.
        auto values = apply_visitor(
            [&](auto const& control) {
                return control(refrigerant.reference_temperature,
                               0. /* initial time */);
            },
            flowAndTemperatureControl);
        updateHeatTransferCoefficients(values.flow_rate);
    }

    static constexpr int number_of_unknowns = 4;

    std::array<double, number_of_unknowns> pipeHeatCapacities() const
    {
        double const& rho_r = refrigerant.density;
        double const& specific_heat_capacity =
            refrigerant.specific_heat_capacity;
        double const& rho_g = grout.rho_g;
        double const& porosity_g = grout.porosity_g;
        double const& heat_cap_g = grout.heat_cap_g;

        return {{/*i1*/ rho_r * specific_heat_capacity,
                 /*o1*/ rho_r * specific_heat_capacity,
                 /*g1*/ (1.0 - porosity_g) * rho_g * heat_cap_g,
                 /*g2*/ (1.0 - porosity_g) * rho_g * heat_cap_g}};
    }

    std::array<double, number_of_unknowns> pipeHeatConductions() const
    {
        double const& lambda_r = refrigerant.thermal_conductivity;
        double const& rho_r = refrigerant.density;
        double const& Cp_r = refrigerant.specific_heat_capacity;
        double const& alpha_L = _pipes.longitudinal_despersion_length;
        double const& porosity_g = grout.porosity_g;
        double const& lambda_g = grout.lambda_g;

        double const velocity_norm = std::abs(_flow_velocity) * std::sqrt(2);

        // Here we calculate the laplace coefficients in the governing
        // equations of BHE. These governing equations can be found in
        // 1) Diersch (2013) FEFLOW book on page 952, M.120-122, or
        // 2) Diersch (2011) Comp & Geosci 37:1122-1135, Eq. 19-22.
        return {{// pipe i1, Eq. 19
                 (lambda_r + rho_r * Cp_r * alpha_L * velocity_norm),
                 // pipe o1, Eq. 20
                 (lambda_r + rho_r * Cp_r * alpha_L * velocity_norm),
                 // pipe g1, Eq. 21
                 (1.0 - porosity_g) * lambda_g,
                 // pipe g2, Eq. 22
                 (1.0 - porosity_g) * lambda_g}};
    }

    std::array<Eigen::Vector3d, number_of_unknowns> pipeAdvectionVectors() const
    {
        double const& rho_r = refrigerant.density;
        double const& Cp_r = refrigerant.specific_heat_capacity;

        return {{// pipe i1, Eq. 19
                 {0, 0, -rho_r * Cp_r * _flow_velocity},
                 // pipe o1, Eq. 20
                 {0, 0, rho_r * Cp_r * _flow_velocity},
                 // grout g1, Eq. 21
                 {0, 0, 0},
                 // grout g2, Eq. 22
                 {0, 0, 0}}};
    }

    template <int NPoints,
              typename SingleUnknownMatrixType,
              typename RMatrixType,
              typename RPiSMatrixType,
              typename RSMatrixType>
    void assembleRMatrices(
        int const idx_bhe_unknowns,
        Eigen::MatrixBase<SingleUnknownMatrixType> const& matBHE_loc_R,
        Eigen::MatrixBase<RMatrixType>& R_matrix,
        Eigen::MatrixBase<RPiSMatrixType>& R_pi_s_matrix,
        Eigen::MatrixBase<RSMatrixType>& R_s_matrix) const
    {
        switch (idx_bhe_unknowns)
        {
            case 0:  // PHI_fig
                R_matrix.block(0, 2 * NPoints, NPoints, NPoints) +=
                    -1.0 * matBHE_loc_R;
                R_matrix.block(2 * NPoints, 0, NPoints, NPoints) +=
                    -1.0 * matBHE_loc_R;

                R_matrix.block(0, 0, NPoints, NPoints) +=
                    1.0 * matBHE_loc_R;  // K_i1
                R_matrix.block(2 * NPoints, 2 * NPoints, NPoints, NPoints) +=
                    1.0 * matBHE_loc_R;  // K_ig
                break;
            case 1:  // PHI_fog
                R_matrix.block(NPoints, 3 * NPoints, NPoints, NPoints) +=
                    -1.0 * matBHE_loc_R;
                R_matrix.block(3 * NPoints, NPoints, NPoints, NPoints) +=
                    -1.0 * matBHE_loc_R;

                R_matrix.block(NPoints, NPoints, NPoints, NPoints) +=
                    1.0 * matBHE_loc_R;  // K_o1
                R_matrix.block(3 * NPoints, 3 * NPoints, NPoints, NPoints) +=
                    1.0 * matBHE_loc_R;  // K_og
                break;
            case 2:  // PHI_gg
                R_matrix.block(2 * NPoints, 3 * NPoints, NPoints, NPoints) +=
                    -1.0 * matBHE_loc_R;
                R_matrix.block(3 * NPoints, 2 * NPoints, NPoints, NPoints) +=
                    -1.0 * matBHE_loc_R;

                R_matrix.block(2 * NPoints, 2 * NPoints, NPoints, NPoints) +=
                    1.0 * matBHE_loc_R;  // K_ig  // notice we only have
                                         // 1 PHI_gg term here.
                R_matrix.block(3 * NPoints, 3 * NPoints, NPoints, NPoints) +=
                    1.0 * matBHE_loc_R;  // K_og  // see Diersch 2013 FEFLOW
                                         // book page 954 Table M.2
                break;
            case 3:  // PHI_gs
                R_s_matrix.template block<NPoints, NPoints>(0, 0).noalias() +=
                    1.0 * matBHE_loc_R;

                R_pi_s_matrix.block(2 * NPoints, 0, NPoints, NPoints) +=
                    -1.0 * matBHE_loc_R;
                R_pi_s_matrix.block(3 * NPoints, 0, NPoints, NPoints) +=
                    -1.0 * matBHE_loc_R;
                R_matrix.block(2 * NPoints, 2 * NPoints, NPoints, NPoints) +=
                    1.0 * matBHE_loc_R;  // K_ig
                R_matrix.block(3 * NPoints, 3 * NPoints, NPoints, NPoints) +=
                    1.0 * matBHE_loc_R;  // K_og
                break;
        }
    }

    /**
     * return the inflow temperature based on outflow temperature and fixed
     * power.
     */
    double getTinByTout(double T_out, double current_time);

    double thermalResistance(int const unknown_index) const
    {
        return _thermal_resistances[unknown_index];
    }

    static constexpr std::pair<int, int> inflow_outflow_bc_component_ids[] = {
        {0, 1}};

private:
    // Placing it here before using it in the cross_section_areas.
    PipeConfiguration1U const _pipes;

public:
    std::array<double, number_of_unknowns> const cross_section_areas = {
        {_pipes.inlet.area(), _pipes.inlet.area(),
         borehole_geometry.area() / 2 - _pipes.outlet.area(),
         borehole_geometry.area() / 2 - _pipes.outlet.area()}};

private:
    void updateHeatTransferCoefficients(double const flow_rate);

    std::array<double, number_of_unknowns> calcThermalResistances(
        double const Nu);

private:
    /// PHI_fig, PHI_fog, PHI_gg, PHI_gs;
    /// Here we store the thermal resistances needed for computation of the heat
    /// exchange coefficients in the governing equations of BHE.
    /// These governing equations can be found in
    /// 1) Diersch (2013) FEFLOW book on page 958, M.3, or
    /// 2) Diersch (2011) Comp & Geosci 37:1122-1135, Eq. 90-97.
    std::array<double, number_of_unknowns> _thermal_resistances;

    /// Flow velocity inside the pipes. Depends on the flow_rate.
    double _flow_velocity;
};
}  // namespace BHE
}  // namespace HeatTransportBHE
}  // namespace ProcessLib
