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

#include "BHECommonCoaxial.h"
#include "FlowAndTemperatureControl.h"
#include "PipeConfigurationCoaxial.h"

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
class BHE_CXC final : public BHECommonCoaxial
{
public:
    BHE_CXC(BoreholeGeometry const& borehole,
            RefrigerantProperties const& refrigerant,
            GroutParameters const& grout,
            FlowAndTemperatureControl const& flowAndTemperatureControl,
            PipeConfigurationCoaxial const& pipes);

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

public:
    std::array<double, number_of_unknowns> const cross_section_areas = {
        {_pipes.inner_pipe.area(),
         _pipes.outer_pipe.area() - _pipes.inner_pipe.outsideArea(),
         borehole_geometry.area() - _pipes.outer_pipe.outsideArea()}};

private:
    std::array<double, number_of_unknowns> calcThermalResistances(
        double const Nu_inner_pipe, double const Nu_annulus_pipe) override;

    std::array<double, 2> velocities() const override
    {
        return {_flow_velocity_inner, _flow_velocity_annulus};
    }
};
}  // namespace BHE
}  // namespace HeatTransportBHE
}  // namespace ProcessLib
