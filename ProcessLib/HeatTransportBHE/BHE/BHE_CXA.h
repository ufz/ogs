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
 * The BHE_CXA class is the realization of Coaxial pipe with Annular type of the
 * Borehole Heate Exchanger. In this class, the pipe heat capacity,
 * pipe heat conduction, pipe advection vectors are intialized according to the
 * geometry of CXA type of BHE. For CXA type of BHE, 3 primary unknowns are
 * assigned on the 1D BHE elements. They are the temperature in inflow pipe
 * T_in, temperature in outflow pipe T_out, temperature of the grout zone
 * surrounding the inflow pipe T_g. These primary variables are solved
 * according to heat convection and conduction equations on the pipes and also
 * in the grout zone. The interaction of the CXA type of BHE and the
 * surrounding soil is regulated through the thermal resistance values, which
 * are calculated specifically during the initialization of the class.
 */

class BHE_CXA final : public BHECommonCoaxial
{
public:
    BHE_CXA(BoreholeGeometry const& borehole,
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
            case 0:  // PHI_fig
                R_matrix.block(0, 2 * NPoints, NPoints, NPoints) +=
                    -1.0 * matBHE_loc_R;
                R_matrix.block(2 * NPoints, 0, NPoints, NPoints) +=
                    -1.0 * matBHE_loc_R;

                R_matrix.block(0, 0, NPoints, NPoints) +=
                    1.0 * matBHE_loc_R;  // K_i
                R_matrix.block(2 * NPoints,
                               2 * NPoints,
                               NPoints,
                               NPoints) += 1.0 * matBHE_loc_R;  // K_ig
                return;
            case 1:  // PHI_ff
                R_matrix.block(0, NPoints, NPoints, NPoints) +=
                    -1.0 * matBHE_loc_R;
                R_matrix.block(NPoints, 0, NPoints, NPoints) +=
                    -1.0 * matBHE_loc_R;

                R_matrix.block(0, 0, NPoints,
                               NPoints) += 1.0 * matBHE_loc_R;  // K_i1
                R_matrix.block(NPoints, NPoints, NPoints, NPoints) +=
                    1.0 * matBHE_loc_R;  // K_o
                return;
            case 2:  // PHI_gs
                R_s_matrix += matBHE_loc_R;

                R_pi_s_matrix.block(2 * NPoints, 0, NPoints, NPoints) +=
                    -1.0 * matBHE_loc_R;

                R_matrix.block(2 * NPoints, 2 * NPoints, NPoints,
                               NPoints) += matBHE_loc_R;  // K_ig
                return;
            default:
                OGS_FATAL(
                    "Error!!! In the function BHE_CXA::assembleRMatrices, "
                    "the index of bhe unknowns is out of range! ");
        }
    }

public:
    std::array<double, number_of_unknowns> const cross_section_areas = {
        {_pipes.outer_pipe.area() - _pipes.inner_pipe.outsideArea(),
         _pipes.inner_pipe.area(),
         borehole_geometry.area() - _pipes.outer_pipe.outsideArea()}};

private:
    std::array<double, number_of_unknowns> calcThermalResistances(
        double const Nu_inner_pipe, double const Nu_annulus_pipe) override;

    std::array<double, 2> velocities() const override
    {
        return {_flow_velocity_annulus, _flow_velocity_inner};
    }
};
}  // namespace BHE
}  // namespace HeatTransportBHE
}  // namespace ProcessLib
