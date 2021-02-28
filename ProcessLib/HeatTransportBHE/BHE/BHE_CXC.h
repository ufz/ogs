/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <Eigen/Eigen>

#include "BHECommonCoaxial.h"

namespace ProcessLib
{
namespace HeatTransportBHE
{
namespace BHE
{
/**
 * The BHE_CXC class is the realization of Coaxial pipe with Centred type of the
 * Borehole Heate Exchanger. In this class, the pipe heat capacity,
 * pipe heat conduction, pipe advection vectors are initialized according to the
 * geometry of CXC type of BHE. For CXC type of BHE, 3 primary unknowns are
 * assigned on the 1D BHE elements. They are the temperature in inflow pipe
 * T_in, temperature in outflow pipe T_out, temperature of the grout zone
 * surrounding the outflow pipe T_g. These primary variables are solved
 * according to heat convection and conduction equations on the pipes and also
 * in the grout zone. The interaction of the CXC type of BHE and the
 * surrounding soil is regulated through the thermal resistance values, which
 * are calculated specifically during the initialization of the class.
 */
class BHE_CXC final : public BHECommonCoaxial
{
public:
    BHE_CXC(BoreholeGeometry const& borehole,
            RefrigerantProperties const& refrigerant,
            GroutParameters const& grout,
            FlowAndTemperatureControl const& flowAndTemperatureControl,
            PipeConfigurationCoaxial const& pipes,
            bool const use_python_bcs)
        : BHECommonCoaxial{borehole, refrigerant,
                           grout,    flowAndTemperatureControl,
                           pipes,    use_python_bcs}
    {
        // Initialize thermal resistances.
        auto values = visit(
            [&](auto const& control) {
                return control(refrigerant.reference_temperature,
                               0. /* initial time */);
            },
            flowAndTemperatureControl);
        updateHeatTransferCoefficients(values.flow_rate);
    }

    template <int NPoints, typename SingleUnknownMatrixType,
              typename RMatrixType, typename RPiSMatrixType,
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

    std::array<double, number_of_unknowns> crossSectionAreas() const
    {
        return {cross_section_area_inner_pipe, cross_section_area_annulus,
                cross_section_area_grout};
    }

private:
    std::array<double, 2> velocities() const override
    {
        return {_flow_velocity_inner, _flow_velocity_annulus};
    }

    std::array<double, number_of_unknowns> getThermalResistances(
        double const& R_gs, double const& R_ff,
        double const& R_fg) const override
    {
        return {R_ff, R_fg, R_gs};
    }
};
}  // namespace BHE
}  // namespace HeatTransportBHE
}  // namespace ProcessLib
