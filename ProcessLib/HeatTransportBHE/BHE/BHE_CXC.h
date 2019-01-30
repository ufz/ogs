/**
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
#include "FlowAndTemperatureControl.h"
#include "PipeConfigurationCXC.h"

namespace ProcessLib
{
namespace HeatTransportBHE
{
namespace BHE
{
/**
 * The BHE_CXC class is the realization of Coaxial pipe with Centred type of the
 * Borehole Heate Exchanger. In this class, the pipe heat capacity,
 * pipe heat conduction, pipe advection vectors are intialized according to the
 * geometry of CXC type of BHE. For CXC type of BHE, 3 primary unknowns are
 * assigned on the 1D BHE elements. They are the temperature in inflow pipe
 * T_in, temperature in outflow pipe T_out, temperature of the grout zone
 * surrounding the outflow pipe T_g. These primary variables are solved
 * according to heat convection and conduction equations on the pipes and also
 * in the grout zone. The interaction of the CXC type of BHE and the
 * surrounding soil is regulated through the thermal resistance values, which
 * are calculated specifically during the initialization of the class.
 */
class BHE_CXC final : public BHECommon
{
public:
    BHE_CXC(BoreholeGeometry const& borehole,
            RefrigerantProperties const& refrigerant,
            GroutParameters const& grout,
            FlowAndTemperatureControl const& flowAndTemperatureControl,
            PipeConfigurationCXC const& pipes);

    static constexpr int number_of_unknowns = 3;
    static constexpr int number_of_grout_zones = 1;

    std::array<double, number_of_unknowns> pipeHeatCapacities() const;

    std::array<double, number_of_unknowns> pipeHeatConductions() const;

    std::array<Eigen::Vector3d, number_of_unknowns> pipeAdvectionVectors()
        const;

    template <int NPoints, typename SingleUnknownMatrixType,
              typename RMatrixType, typename RPiSMatrixType,
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
            case 0:  // PHI_ff
                R_matrix.block(0, NPoints, NPoints, NPoints) +=
                    -1.0 * matBHE_loc_R;
                R_matrix.block(NPoints, 0, NPoints, NPoints) +=
                    -1.0 * matBHE_loc_R;

                R_matrix.block(0, 0, NPoints,
                               NPoints) += 1.0 * matBHE_loc_R;  // K_i
                R_matrix.block(NPoints, NPoints, NPoints, NPoints) +=
                    1.0 * matBHE_loc_R;  // K_o
                return;
            case 1:  // PHI_fog
                R_matrix.block(NPoints, 2 * NPoints, NPoints, NPoints) +=
                    -1.0 * matBHE_loc_R;
                R_matrix.block(2 * NPoints, NPoints, NPoints, NPoints) +=
                    -1.0 * matBHE_loc_R;

                R_matrix.block(NPoints, NPoints, NPoints, NPoints) +=
                    1.0 * matBHE_loc_R;  // K_o
                R_matrix.block(2 * NPoints,
                               2 * NPoints,
                               NPoints,
                               NPoints) += 1.0 * matBHE_loc_R;  // K_og
                return;
            case 2:  // PHI_gs
                R_s_matrix += matBHE_loc_R;

                R_pi_s_matrix.block(2 * NPoints, 0, NPoints, NPoints) +=
                    -1.0 * matBHE_loc_R;

                R_matrix.block(2 * NPoints, 2 * NPoints, NPoints,
                               NPoints) += matBHE_loc_R;  // K_og
                return;
            default:
                OGS_FATAL(
                    "Error!!! In the function BHE_CXC::assembleRMatrices, "
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

private:
    // Placing it here before using it in the cross_section_areas.
    PipeConfigurationCXC const _pipes;

public:
    std::array<double, number_of_unknowns> const cross_section_areas = {
        {_pipes.inner_inflow_pipe.area(),
         _pipes.outer_pipe.area() - _pipes.inner_inflow_pipe.outsideArea(),
         borehole_geometry.area() - _pipes.outer_pipe.outsideArea()}};

private:
    void updateHeatTransferCoefficients(double const flow_rate);

    std::array<double, number_of_unknowns> calcThermalResistances(
        double const Nu_o, double const Nu_i);

private:
    /// PHI_ff, PHI_fog, PHI_gs;
    /// Here we store the thermal resistances needed for computation of the heat
    /// exchange coefficients in the governing equations of BHE.
    /// These governing equations can be found in
    /// 1) Diersch (2013) FEFLOW book on page 958, M.3, or
    /// 2) Diersch (2011) Comp & Geosci 37:1122-1135, Eq. 90-97.
    std::array<double, number_of_unknowns> _thermal_resistances;

    /// Flow velocity inside the pipes and annulus. Depends on the flow_rate.
    double _flow_velocity, _flow_velocity_annulus;
};
}  // namespace BHE
}  // namespace HeatTransportBHE
}  // namespace ProcessLib
