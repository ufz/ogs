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

#include <Eigen/Eigen>

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
 * The BHE_1P class is the realization of single-pipe type of Borehole Heate
 * Exchanger. In this class, the pipe heat capacity, pipe heat conductiion, pie
 * advection vectors are intialized according to the geometry of the single-pipe
 * type of BHE. For this type of BHE, 2 primary unknowns are assigned on the 1D
 * BHE elements. They are the temperature in the pipe T_p, and temperature of
 * the grout zone sorrounding the single pipe T_g. These two primary varaibles
 * are solved according to heat convection and conduction equations on the pipes
 * and also in the grout zones. The interaction of the 1P type of BHE and the
 * sorrounding soil is regulated through the thermal resistance values, which
 * are calculated specifically during the initialization of the class.
 */
class BHE_1P final : public BHECommonUType
{
public:
    BHE_1P(BoreholeGeometry const& borehole,
           RefrigerantProperties const& refrigerant,
           GroutParameters const& grout,
           FlowAndTemperatureControl const& flowAndTemperatureControl,
           PipeConfigurationUType const& pipes)
        : BHECommonUType{borehole, refrigerant, grout,
                         flowAndTemperatureControl, pipes}
    {
        _thermal_resistances.fill(std::numeric_limits<double>::quiet_NaN());

        // Initialize thermal resistances.
        auto values = visit(
            [&](auto const& control) {
                return control(refrigerant.reference_temperature,
                               0. /* initial time */);
            },
            flowAndTemperatureControl);
        updateHeatTransferCoefficients(values.flow_rate);
    }

    static constexpr int number_of_unknowns = 2;
    static constexpr int number_of_grout_zones = 1;

    std::array<double, number_of_unknowns> pipeHeatCapacities() const;

    std::array<double, number_of_unknowns> pipeHeatConductions() const;

    std::array<Eigen::Vector3d, number_of_unknowns> pipeAdvectionVectors()
        const;

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
        // TODO, this needs to be changed.
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
                return;
            case 1:  // PHI_fog
                R_matrix.block(NPoints, 3 * NPoints, NPoints, NPoints) +=
                    -1.0 * matBHE_loc_R;
                R_matrix.block(3 * NPoints, NPoints, NPoints, NPoints) +=
                    -1.0 * matBHE_loc_R;

                R_matrix.block(NPoints, NPoints, NPoints, NPoints) +=
                    1.0 * matBHE_loc_R;  // K_o1
                R_matrix.block(3 * NPoints, 3 * NPoints, NPoints, NPoints) +=
                    1.0 * matBHE_loc_R;  // K_og
                return;
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
                return;
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
                return;
            default:
                OGS_FATAL(
                    "Error!!! In the function BHE_1U::assembleRMatrices, "
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
        {0, 1}};

public:
    std::array<double, number_of_unknowns> crossSectionAreas() const
    {
        return {{_pipes.inlet.area(),
                 borehole_geometry.area() - _pipes.inlet.area()}};
    }

private:
    void updateHeatTransferCoefficients(double const flow_rate);

    std::array<double, number_of_unknowns> calcThermalResistances(
        double const Nu);

    double compute_R_gs(double const chi, double const R_g);

private:
    /// PHI_fig, PHI_gs;
    /// Here we store the thermal resistances needed for computation of the heat
    /// exchange coefficients in the governing equations of BHE.
    /// These governing equations can be found in
    /// 1) Diersch (2013) FEFLOW book on page 958, M.3, or
    /// 2) Diersch (2011) Comp & Geosci 37:1122-1135, Eq. 90-97.
    std::array<double, number_of_unknowns> _thermal_resistances;
};
}  // namespace BHE
}  // namespace HeatTransportBHE
}  // namespace ProcessLib
