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
#include <optional>

#include "BaseLib/Error.h"

#include "BHECommon.h"
#include "FlowAndTemperatureControl.h"
#include "PipeConfiguration1PType.h"

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
class BHE_1P final : public BHECommon
{
public:
    BHE_1P(BoreholeGeometry const& borehole,
           RefrigerantProperties const& refrigerant,
           GroutParameters const& grout,
           FlowAndTemperatureControl const& flowAndTemperatureControl,
           PipeConfiguration1PType const& pipes,
           bool const use_python_bcs);

    static constexpr int number_of_unknowns = 2;
    static constexpr int number_of_grout_zones = 1;

    std::array<double, number_of_unknowns> pipeHeatCapacities() const;

    std::array<double, number_of_unknowns> pipeHeatConductions() const;

    std::array<Eigen::Vector3d, number_of_unknowns> pipeAdvectionVectors(
        Eigen::Vector3d const& elem_direction) const;

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
        // Here we are looping over two resistance terms
        // First PHI_fg is the resistance between pipe and grout
        // Second PHI_gs is the resistance between grout and soil
        switch (idx_bhe_unknowns)
        {
            case 0:  // PHI_fg
                R_matrix.block(0, NPoints, NPoints, NPoints) +=
                    -1.0 * matBHE_loc_R;
                R_matrix.block(NPoints, 0, NPoints, NPoints) +=
                    -1.0 * matBHE_loc_R;

                R_matrix.block(0, 0, NPoints, NPoints) +=
                    matBHE_loc_R;  // K_i/o
                R_matrix.block(NPoints, NPoints, NPoints, NPoints) +=
                    matBHE_loc_R;  // K_fg
                return;
            case 1:  // PHI_gs
                R_s_matrix += matBHE_loc_R;

                R_pi_s_matrix.block(NPoints, 0, NPoints, NPoints) +=
                    -1.0 * matBHE_loc_R;

                R_matrix.block(NPoints, NPoints, NPoints, NPoints) +=
                    matBHE_loc_R;  // K_fg
                return;
            default:
                OGS_FATAL(
                    "Error!!! In the function BHE_1P::assembleRMatrices, "
                    "the index of bhe resistance term is out of range! ");
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

    static std::array<std::pair<std::size_t /*node_id*/, int /*component*/>, 2>
    getBHEInflowDirichletBCNodesAndComponents(std::size_t const top_node_id,
                                              std::size_t const bottom_node_id,
                                              int const in_component_id);

    static std::optional<
        std::array<std::pair<std::size_t /*node_id*/, int /*component*/>, 2>>
    getBHEBottomDirichletBCNodesAndComponents(
        std::size_t const /*bottom_node_id*/,
        int const /*in_component_id*/,
        int const /*out_component_id*/);

public:
    std::array<double, number_of_unknowns> crossSectionAreas() const;

    void updateHeatTransferCoefficients(double const flow_rate);

protected:
    PipeConfiguration1PType const _pipe;

    /// Flow velocity inside the pipes. Depends on the flow_rate.
    double _flow_velocity = std::numeric_limits<double>::quiet_NaN();

private:
    std::array<double, number_of_unknowns> calcThermalResistances(
        double const Nu);

    static double compute_R_gs(double const chi, double const R_g);

private:
    /// PHI_fg, PHI_gs;
    /// Here we store the thermal resistances needed for computation of the heat
    /// exchange coefficients in the governing equations of BHE.
    std::array<double, number_of_unknowns> _thermal_resistances;
};
}  // namespace BHE
}  // namespace HeatTransportBHE
}  // namespace ProcessLib
