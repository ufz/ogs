// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <Eigen/Core>

#include "BHECommonCoaxial.h"
#include "BHEParameterValidation.h"
#include "BaseLib/Error.h"
#include "ParameterLib/SpatialPosition.h"

namespace ProcessLib
{
namespace HeatTransportBHE
{
namespace BHE
{
class BHE_CXA final : public BHECommonCoaxial
{
public:
    BHE_CXA(BoreholeGeometry const& borehole,
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
            [&](auto const& control)
            {
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
                    "BHE_CXA::assembleRMatrices: unknown index {:d} "
                    "out of range.",
                    idx_bhe_unknowns);
        }
    }

    std::array<double, number_of_unknowns> crossSectionAreas(
        ParameterLib::SpatialPosition const& pos) const
    {
        double const D = sampleStrictPositive(borehole_geometry.diameter, 0.0,
                                              pos, "borehole_diameter");
        double const borehole_area = Pipe::circleArea(D);
        return {cross_section_area_annulus, cross_section_area_inner_pipe,
                checkedGroutArea(borehole_area, _pipes.outer_pipe.outsideArea(),
                                 pos)};
    }

private:
    void assignVelocities(double inner_vel, double annulus_vel) override
    {
        // CXA: unknown 0 = annulus (inflow), unknown 1 = inner pipe (outflow).
        // velocity_inner stores the velocity for unknown 0, so it must
        // receive the annulus velocity here.
        velocity_inner = annulus_vel;
        velocity_annulus = inner_vel;
    }

    std::vector<double> getThermalResistances(double const& R_gs,
                                              double const& R_ff,
                                              double const& R_fg) const override
    {
        return {R_fg, R_ff, R_gs};
    }
};
}  // namespace BHE
}  // namespace HeatTransportBHE
}  // namespace ProcessLib
