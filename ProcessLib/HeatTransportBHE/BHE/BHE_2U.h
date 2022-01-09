/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <Eigen/Eigen>

#include <optional>
#include "BaseLib/Error.h"

#include "BHECommon.h"
#include "BHECommonUType.h"
#include "FlowAndTemperatureControl.h"
#include "PipeConfigurationUType.h"

namespace ProcessLib
{
namespace HeatTransportBHE
{
namespace BHE
{
/**
 * The BHE_2U class is the realization of 2U type of Borehole Heate Exchanger.
 * In this class, the pipe heat capacity, pipe heat conductiion, pipe advection
 * vectors are initialized according to the geometry of 2U type of BHE.
 * For 2U type of BHE, 8 primary unknowns are assigned on the 1D BHE elements.
 * They are the temperature in inflow pipe T_i1 and T_i2, temperature in outflow
 * pipe T_o1 and T_o2, temperature of the four grout zones surrounding the
 * inflow and outflow pipe T_g1, T_g2, T_g3 and T_g4. These primary variables
 * are solved according to heat convection and conduction equations on the pipes
 * and also in the grout zones. The interaction of the 2U type of BHE and the
 * surrounding soil is regulated through the thermal resistance values, which
 * are calculated specifically during the initialization of the class.
 */
class BHE_2U final : public BHECommonUType
{
public:
    BHE_2U(BoreholeGeometry const& borehole,
           RefrigerantProperties const& refrigerant,
           GroutParameters const& grout,
           FlowAndTemperatureControl const& flowAndTemperatureControl,
           PipeConfigurationUType const& pipes,
           bool const use_python_bcs);

    static constexpr int number_of_unknowns = 8;
    static constexpr int number_of_grout_zones = 4;

    std::array<double, number_of_unknowns> pipeHeatCapacities() const;

    std::array<double, number_of_unknowns> pipeHeatConductions() const;

    std::array<Eigen::Vector3d, number_of_unknowns> pipeAdvectionVectors(
        Eigen::Vector3d const& /*elem_direction*/) const;

    template <int NPoints,
              typename SingleUnknownMatrixType,
              typename RMatrixType,
              typename RPiSMatrixType,
              typename RSMatrixType>
    static void assembleRMatrices(
        int const idx_bhe_unknowns,
        Eigen::MatrixBase<SingleUnknownMatrixType> const& matBHE_loc_R,
        Eigen::MatrixBase<RMatrixType>& R_matrix,
        Eigen::MatrixBase<RPiSMatrixType>& R_pi_s_matrix,
        Eigen::MatrixBase<RSMatrixType>& R_s_matrix)
    {
        switch (idx_bhe_unknowns)
        {
            case 0:  // PHI_fig
                R_matrix.block(0, 4 * NPoints, NPoints, NPoints) +=
                    -1.0 * matBHE_loc_R;
                R_matrix.block(4 * NPoints, 0, NPoints, NPoints) +=
                    -1.0 * matBHE_loc_R;  // R_i1

                R_matrix.block(2, 5 * NPoints, NPoints, NPoints) +=
                    -1.0 * matBHE_loc_R;
                R_matrix.block(5 * NPoints, 2, NPoints, NPoints) +=
                    -1.0 * matBHE_loc_R;  // R_i2

                R_matrix.block(0, 0, NPoints, NPoints) +=
                    1.0 * matBHE_loc_R;  // K_i1
                R_matrix.block(NPoints, NPoints, NPoints, NPoints) +=
                    1.0 * matBHE_loc_R;  // K_i2
                R_matrix.block(4 * NPoints, 4 * NPoints, NPoints, NPoints) +=
                    1.0 * matBHE_loc_R;
                R_matrix.block(5 * NPoints, 5 * NPoints, NPoints, NPoints) +=
                    1.0 * matBHE_loc_R;  // K_ig
                return;
            case 1:  // PHI_fog
                R_matrix.block(2 * NPoints, 6 * NPoints, NPoints, NPoints) +=
                    -1.0 * matBHE_loc_R;
                R_matrix.block(6 * NPoints, 2 * NPoints, NPoints, NPoints) +=
                    -1.0 * matBHE_loc_R;  // R_o1
                R_matrix.block(3 * NPoints, 7 * NPoints, NPoints, NPoints) +=
                    -1.0 * matBHE_loc_R;
                R_matrix.block(7 * NPoints, 3 * NPoints, NPoints, NPoints) +=
                    -1.0 * matBHE_loc_R;  // R_o2

                R_matrix.block(2 * NPoints, 2 * NPoints, NPoints, NPoints) +=
                    1.0 * matBHE_loc_R;  // K_o1
                R_matrix.block(3 * NPoints, 3 * NPoints, NPoints, NPoints) +=
                    1.0 * matBHE_loc_R;  // K_o2
                R_matrix.block(6 * NPoints, 6 * NPoints, NPoints, NPoints) +=
                    1.0 * matBHE_loc_R;
                R_matrix.block(7 * NPoints, 7 * NPoints, NPoints, NPoints) +=
                    1.0 * matBHE_loc_R;  // K_og
                return;
            case 2:  // PHI_gg_1
                R_matrix.block(4 * NPoints, 6 * NPoints, NPoints, NPoints) +=
                    -1.0 * matBHE_loc_R;
                R_matrix.block(6 * NPoints, 4 * NPoints, NPoints, NPoints) +=
                    -1.0 * matBHE_loc_R;
                R_matrix.block(4 * NPoints, 7 * NPoints, NPoints, NPoints) +=
                    -1.0 * matBHE_loc_R;
                R_matrix.block(7 * NPoints, 4 * NPoints, NPoints, NPoints) +=
                    -1.0 * matBHE_loc_R;
                R_matrix.block(5 * NPoints, 6 * NPoints, NPoints, NPoints) +=
                    -1.0 * matBHE_loc_R;
                R_matrix.block(6 * NPoints, 5 * NPoints, NPoints, NPoints) +=
                    -1.0 * matBHE_loc_R;
                R_matrix.block(5 * NPoints, 7 * NPoints, NPoints, NPoints) +=
                    -1.0 * matBHE_loc_R;
                R_matrix.block(7 * NPoints, 5 * NPoints, NPoints, NPoints) +=
                    -1.0 * matBHE_loc_R;  // R_g1

                R_matrix.block(4 * NPoints, 4 * NPoints, NPoints, NPoints) +=
                    2.0 * matBHE_loc_R;
                R_matrix.block(5 * NPoints, 5 * NPoints, NPoints, NPoints) +=
                    2.0 * matBHE_loc_R;  // K_ig  // notice we have two
                                         // PHI_gg_1 term here.
                R_matrix.block(6 * NPoints, 6 * NPoints, NPoints, NPoints) +=
                    2.0 * matBHE_loc_R;
                R_matrix.block(7 * NPoints, 7 * NPoints, NPoints, NPoints) +=
                    2.0 * matBHE_loc_R;  // K_og  // see Diersch 2013 FEFLOW
                                         // book page 954 Table M.2
                return;
            case 3:  // PHI_gg_2
                R_matrix.block(4 * NPoints, 5 * NPoints, NPoints, NPoints) +=
                    -1.0 * matBHE_loc_R;
                R_matrix.block(5 * NPoints, 4 * NPoints, NPoints, NPoints) +=
                    -1.0 * matBHE_loc_R;
                R_matrix.block(6 * NPoints, 7 * NPoints, NPoints, NPoints) +=
                    -1.0 * matBHE_loc_R;
                R_matrix.block(7 * NPoints, 6 * NPoints, NPoints, NPoints) +=
                    -1.0 * matBHE_loc_R;  // R_g2

                R_matrix.block(4 * NPoints, 4 * NPoints, NPoints, NPoints) +=
                    1.0 * matBHE_loc_R;
                R_matrix.block(5 * NPoints, 5 * NPoints, NPoints, NPoints) +=
                    1.0 * matBHE_loc_R;  // K_ig  // notice we only have
                                         // 1 PHI_gg term here.
                R_matrix.block(6 * NPoints, 6 * NPoints, NPoints, NPoints) +=
                    1.0 * matBHE_loc_R;
                R_matrix.block(7 * NPoints, 7 * NPoints, NPoints, NPoints) +=
                    1.0 * matBHE_loc_R;  // K_og  // see Diersch 2013 FEFLOW
                                         // book page 954 Table M.2
                return;
            case 4:  // PHI_gs
                R_s_matrix.template block<NPoints, NPoints>(0, 0).noalias() +=
                    1.0 * matBHE_loc_R;

                R_pi_s_matrix.block(4 * NPoints, 0, NPoints, NPoints) +=
                    -1.0 * matBHE_loc_R;
                R_pi_s_matrix.block(5 * NPoints, 0, NPoints, NPoints) +=
                    -1.0 * matBHE_loc_R;
                R_pi_s_matrix.block(6 * NPoints, 0, NPoints, NPoints) +=
                    -1.0 * matBHE_loc_R;
                R_pi_s_matrix.block(7 * NPoints, 0, NPoints, NPoints) +=
                    -1.0 * matBHE_loc_R;
                R_matrix.block(4 * NPoints, 4 * NPoints, NPoints, NPoints) +=
                    1.0 * matBHE_loc_R;
                R_matrix.block(5 * NPoints, 5 * NPoints, NPoints, NPoints) +=
                    1.0 * matBHE_loc_R;  // K_ig
                R_matrix.block(6 * NPoints, 6 * NPoints, NPoints, NPoints) +=
                    1.0 * matBHE_loc_R;
                R_matrix.block(7 * NPoints, 7 * NPoints, NPoints, NPoints) +=
                    1.0 * matBHE_loc_R;  // K_og
                return;
            default:
                OGS_FATAL(
                    "Error!!! In the function BHE_2U::assembleRMatrices, "
                    "the index of bhe unknowns is out of range! ");
        }
    }

    /// Return the inflow temperature for the boundary condition.
    double updateFlowRateAndTemperature(double T_out, double current_time);

    double thermalResistance(int const unknown_index) const
    {
        return _thermal_resistances[unknown_index];
    }

    static constexpr std::pair<int, int> inflow_outflow_bc_component_ids[] = {
        {0, 2}, {1, 3}};

    static std::array<std::pair<std::size_t /*node_id*/, int /*component*/>, 2>
    getBHEInflowDirichletBCNodesAndComponents(
        std::size_t const top_node_id,
        std::size_t const /*bottom_node_id*/,
        int const in_component_id);

    static std::optional<
        std::array<std::pair<std::size_t /*node_id*/, int /*component*/>, 2>>
    getBHEBottomDirichletBCNodesAndComponents(std::size_t const bottom_node_id,
                                              int const in_component_id,
                                              int const out_component_id);

public:
    std::array<double, number_of_unknowns> crossSectionAreas() const;

    void updateHeatTransferCoefficients(double const flow_rate);

private:
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
