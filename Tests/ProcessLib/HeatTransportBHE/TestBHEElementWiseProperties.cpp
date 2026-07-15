// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <Eigen/Core>
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <numbers>
#include <stdexcept>
#include <vector>

#include "BaseLib/ConfigTree.h"
#include "MeshLib/Elements/Line.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"
#include "MeshLib/Utils/addPropertyToMesh.h"
#include "ParameterLib/ConstantParameter.h"
#include "ParameterLib/SpatialPosition.h"
#include "ProcessLib/HeatTransportBHE/BHE/BHE_1P.h"
#include "ProcessLib/HeatTransportBHE/BHE/BHE_1U.h"
#include "ProcessLib/HeatTransportBHE/BHE/BHE_2U.h"
#include "ProcessLib/HeatTransportBHE/BHE/BHE_CXA.h"
#include "ProcessLib/HeatTransportBHE/BHE/BHE_CXC.h"
#include "ProcessLib/HeatTransportBHE/BHE/BoreholeGeometry.h"
#include "ProcessLib/HeatTransportBHE/BHE/FlowAndTemperatureControl.h"
#include "ProcessLib/HeatTransportBHE/BHE/GroutParameters.h"
#include "ProcessLib/HeatTransportBHE/BHE/MeshUtils.h"
#include "ProcessLib/HeatTransportBHE/BHE/Pipe.h"
#include "ProcessLib/HeatTransportBHE/BHE/PipeConfiguration1PType.h"
#include "ProcessLib/HeatTransportBHE/BHE/PipeConfigurationCoaxial.h"
#include "ProcessLib/HeatTransportBHE/BHE/PipeConfigurationUType.h"
#include "ProcessLib/HeatTransportBHE/BHE/RefrigerantProperties.h"
#include "Tests/TestTools.h"

namespace
{
namespace BHEns = ProcessLib::HeatTransportBHE::BHE;

/// Helper to build a SpatialPosition for a given element ID.
ParameterLib::SpatialPosition makePos(std::size_t const element_id)
{
    ParameterLib::SpatialPosition pos;
    pos.setElementID(element_id);
    return pos;
}

/// Parse a small XML snippet into a ConfigTree for factory-level tests.
BaseLib::ConfigTree makeConfigTree(char const* xml)
{
    return BaseLib::ConfigTree(Tests::readXml(xml),
                               "TestBHEElementWiseProperties",
                               BaseLib::ConfigTree::onerror,
                               BaseLib::ConfigTree::onwarning);
}

// RefrigerantProperties fields: dynamic_viscosity [kg m^-1 s^-1],
// density [kg m^-3], thermal_conductivity [W m^-1 K^-1],
// specific_heat_capacity [J kg^-1 K^-1], reference_temperature [degC].
BHEns::RefrigerantProperties const kRefrigerant{0.00067418, 992.92, 0.62863,
                                                4198.0, 25.0};
// GroutParameters fields: rho_g (density) [kg m^-3], porosity_g [-],
// heat_cap_g (specific heat capacity) [J kg^-1 K^-1],
// lambda_g (thermal conductivity) [W m^-1 K^-1].
BHEns::GroutParameters const kGrout{2190.0, 0.0, 1735.16, 0.806};

// Standard "valid" parameter values reused across tests. Validation tests
// that need a *bad* diameter or wall TC declare their own local override.
ParameterLib::ConstantParameter<double> const kOkDiameter{"diameter", 0.22};
ParameterLib::ConstantParameter<double> const kOkWallTc{"wall_tc", 0.42};
ParameterLib::ConstantParameter<double> const kFlowRateParam{"flow_rate",
                                                             2.0e-4};
ParameterLib::ConstantParameter<double> const kTemperatureParam{"temperature",
                                                                25.0};

// Geometry of the single-diameter U-type / 1P test pipe. kPipeOutsideDiameter
// is the pipe-center distance at which the R_ar acosh argument is exactly 1
// (pipes touching), so the touching-pipe tests derive their distance from it
// instead of restating the arithmetic.
double const kPipeInsideDiameter = 0.04;
double const kPipeWallThickness = 0.0029;
double const kPipeOutsideDiameter =
    kPipeInsideDiameter + 2 * kPipeWallThickness;  // 0.0458

BHEns::FlowAndTemperatureControl makeConstantControl(
    ParameterLib::Parameter<double> const& temperature_param,
    ParameterLib::Parameter<double> const& flow_rate_param)
{
    // Third argument is flow_rate_min: |flow_rate| below this threshold is
    // clamped to zero. 0.0 means the prescribed flow rate is always used.
    return BHEns::InflowTemperature{temperature_param, flow_rate_param, 0.0};
}
}  // namespace

/// CXA BHE with constant parameters must produce finite thermal resistances.
TEST(BHECommonCoaxial, CxaConstantParamsGiveFiniteResistancesAndPositiveAreas)
{
    using namespace ProcessLib::HeatTransportBHE::BHE;

    BoreholeGeometry const borehole{8.0, kOkDiameter};

    ParameterLib::ConstantParameter<double> const inner_tc{"inner_tc", 0.001};
    ParameterLib::ConstantParameter<double> const outer_tc{"outer_tc", 1.3};
    Pipe const inner_pipe{0.09532, 0.00734, inner_tc};
    Pipe const outer_pipe{0.16626, 0.00587, outer_tc};
    PipeConfigurationCoaxial const pipes{inner_pipe, outer_pipe, 0.001};

    auto const control = makeConstantControl(kTemperatureParam, kFlowRateParam);

    BHE_CXA const bhe{borehole, kRefrigerant, kGrout, control, pipes, false};

    auto const pos = makePos(0);

    // All 3 thermal resistances must be finite.
    for (int u = 0; u < 3; ++u)
    {
        EXPECT_TRUE(std::isfinite(bhe.thermalResistances(pos)[u]))
            << "unknown=" << u;
    }

    // Cross-section areas must all be positive.
    auto const areas = bhe.crossSectionAreas(pos);
    for (int i = 0; i < 3; ++i)
    {
        EXPECT_GT(areas[i], 0.0) << "cross_section_area[" << i << "]";
    }
}

/// BHE_1U with constant parameters must produce finite thermal resistances.
TEST(BHECommonUType, Bhe1UConstantParamsGiveFiniteResistancesAndPositiveAreas)
{
    using namespace ProcessLib::HeatTransportBHE::BHE;

    BoreholeGeometry const borehole{8.0, kOkDiameter};

    Pipe const inlet{kPipeInsideDiameter, kPipeWallThickness, kOkWallTc};
    Pipe const outlet{kPipeInsideDiameter, kPipeWallThickness, kOkWallTc};
    PipeConfigurationUType const pipes{inlet, outlet, 0.06, 0.001};

    auto const control = makeConstantControl(kTemperatureParam, kFlowRateParam);

    BHE_1U const bhe{borehole, kRefrigerant, kGrout, control, pipes, false};

    auto const pos = makePos(0);

    // All 4 thermal resistance unknowns must be finite.
    for (int u = 0; u < 4; ++u)
    {
        EXPECT_TRUE(std::isfinite(bhe.thermalResistances(pos)[u]))
            << "unknown=" << u;
    }

    // Cross-section areas: pipe (0,1) + 2 grout zones (2,3).
    auto const areas = bhe.crossSectionAreas(pos);
    for (int i = 0; i < 4; ++i)
    {
        EXPECT_GT(areas[i], 0.0) << "cross_section_area[" << i << "]";
    }
}

/// BHE_1P with constant parameters must produce finite thermal resistances.
TEST(BHECommon1P, Bhe1PConstantParamsGiveFiniteResistancesAndPositiveAreas)
{
    using namespace ProcessLib::HeatTransportBHE::BHE;

    BoreholeGeometry const borehole{8.0, kOkDiameter};

    Pipe const single_pipe{kPipeInsideDiameter, kPipeWallThickness, kOkWallTc};
    PipeConfiguration1PType const pipes{single_pipe, 0.001};

    auto const control = makeConstantControl(kTemperatureParam, kFlowRateParam);

    BHE_1P const bhe{borehole, kRefrigerant, kGrout, control, pipes, false};

    auto const pos = makePos(0);

    // Both thermal resistance unknowns (pipe-grout, grout-soil) must be finite.
    for (int u = 0; u < 2; ++u)
    {
        EXPECT_TRUE(std::isfinite(bhe.thermalResistances(pos)[u]))
            << "unknown=" << u;
    }

    auto const areas = bhe.crossSectionAreas(pos);
    for (int i = 0; i < 2; ++i)
    {
        EXPECT_GT(areas[i], 0.0) << "cross_section_area[" << i << "]";
    }
}

/// Verify that the ConstantParameter used for borehole diameter produces the
/// expected cross-section area independent of element ID.
TEST(BHECommonCoaxial, ConstantDiameterSameAtAllElements)
{
    using namespace ProcessLib::HeatTransportBHE::BHE;

    BoreholeGeometry const borehole{8.0, kOkDiameter};

    ParameterLib::ConstantParameter<double> const inner_tc{"inner_tc", 0.001};
    ParameterLib::ConstantParameter<double> const outer_tc{"outer_tc", 1.3};
    Pipe const inner_pipe{0.09532, 0.00734, inner_tc};
    Pipe const outer_pipe{0.16626, 0.00587, outer_tc};
    PipeConfigurationCoaxial const pipes{inner_pipe, outer_pipe, 0.001};

    auto const control = makeConstantControl(kTemperatureParam, kFlowRateParam);

    BHE_CXA const bhe{borehole, kRefrigerant, kGrout, control, pipes, false};

    // With a constant diameter, areas and resistances should be the same
    // regardless of element ID.
    auto const pos0 = makePos(0);
    auto const pos5 = makePos(5);

    auto const areas0 = bhe.crossSectionAreas(pos0);
    auto const areas5 = bhe.crossSectionAreas(pos5);
    for (int i = 0; i < 3; ++i)
    {
        EXPECT_DOUBLE_EQ(areas0[i], areas5[i]);
    }

    auto const R0 = bhe.thermalResistances(pos0);
    auto const R5 = bhe.thermalResistances(pos5);
    for (int u = 0; u < 3; ++u)
    {
        EXPECT_DOUBLE_EQ(R0[u], R5[u]);
    }
}

namespace
{
/// Build a BHE_1P with the given diameter and wall thermal conductivity
/// parameters, keeping everything else at standard benchmark values.
BHEns::BHE_1P make1PWith(ParameterLib::Parameter<double> const& diameter_param,
                         ParameterLib::Parameter<double> const& wall_tc_param)
{
    BHEns::BoreholeGeometry const borehole{8.0, diameter_param};
    BHEns::Pipe const single_pipe{kPipeInsideDiameter, kPipeWallThickness,
                                  wall_tc_param};
    BHEns::PipeConfiguration1PType const pipes{single_pipe, 0.001};

    return BHEns::BHE_1P{
        borehole, kRefrigerant,
        kGrout,   makeConstantControl(kTemperatureParam, kFlowRateParam),
        pipes,    false};
}

/// Build a BHE_1U with separate inlet and outlet wall thermal conductivities.
BHEns::BHE_1U make1UWithWallTCs(
    ParameterLib::Parameter<double> const& diameter_param,
    ParameterLib::Parameter<double> const& inlet_wall_tc,
    ParameterLib::Parameter<double> const& outlet_wall_tc,
    double const distance)
{
    BHEns::BoreholeGeometry const borehole{8.0, diameter_param};
    BHEns::Pipe const inlet{kPipeInsideDiameter, kPipeWallThickness,
                            inlet_wall_tc};
    BHEns::Pipe const outlet{kPipeInsideDiameter, kPipeWallThickness,
                             outlet_wall_tc};
    BHEns::PipeConfigurationUType const pipes{inlet, outlet, distance, 0.001};

    return BHEns::BHE_1U{
        borehole, kRefrigerant,
        kGrout,   makeConstantControl(kTemperatureParam, kFlowRateParam),
        pipes,    false};
}

/// Build a BHE_2U with separate inlet and outlet wall thermal conductivities.
BHEns::BHE_2U make2UWithWallTCs(
    ParameterLib::Parameter<double> const& diameter_param,
    ParameterLib::Parameter<double> const& inlet_wall_tc,
    ParameterLib::Parameter<double> const& outlet_wall_tc,
    double const distance)
{
    BHEns::BoreholeGeometry const borehole{8.0, diameter_param};
    BHEns::Pipe const inlet{kPipeInsideDiameter, kPipeWallThickness,
                            inlet_wall_tc};
    BHEns::Pipe const outlet{kPipeInsideDiameter, kPipeWallThickness,
                             outlet_wall_tc};
    BHEns::PipeConfigurationUType const pipes{inlet, outlet, distance, 0.001};

    return BHEns::BHE_2U{
        borehole, kRefrigerant,
        kGrout,   makeConstantControl(kTemperatureParam, kFlowRateParam),
        pipes,    false};
}

/// Build a BHE_1U with the given diameter, wall TC, and pipe-center distance.
BHEns::BHE_1U make1UWith(ParameterLib::Parameter<double> const& diameter_param,
                         ParameterLib::Parameter<double> const& wall_tc_param,
                         double const distance)
{
    return make1UWithWallTCs(diameter_param, wall_tc_param, wall_tc_param,
                             distance);
}

/// Build a BHE_2U with the given diameter, wall TC, and pipe-center distance.
BHEns::BHE_2U make2UWith(ParameterLib::Parameter<double> const& diameter_param,
                         ParameterLib::Parameter<double> const& wall_tc_param,
                         double const distance)
{
    return make2UWithWallTCs(diameter_param, wall_tc_param, wall_tc_param,
                             distance);
}
}  // namespace

/// Factory-level rejection: createPipe() still rejects a non-positive pipe
/// diameter at construction time. This is separate from the wall thermal
/// conductivity, which is now a ParameterLib parameter sampled (and validated
/// for positivity) per element rather than checked entirely at factory time.
TEST(ProcessLibBHEPipe, CreatePipeRejectsInvalidDiameter)
{
    auto cfg = makeConfigTree(
        "<p>"
        "  <diameter>-0.04</diameter>"
        "  <wall_thickness>0.0029</wall_thickness>"
        "  <wall_thermal_conductivity>0.42</wall_thermal_conductivity>"
        "</p>");
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> parameters;
    EXPECT_ANY_THROW(BHEns::createPipe(cfg.getConfigSubtree("p"), parameters));
}

/// Factory-level rejection: createPipe() still rejects a negative pipe wall
/// thickness at construction time.
TEST(ProcessLibBHEPipe, CreatePipeRejectsNegativeWallThickness)
{
    auto cfg = makeConfigTree(
        "<p>"
        "  <diameter>0.04</diameter>"
        "  <wall_thickness>-0.001</wall_thickness>"
        "  <wall_thermal_conductivity>0.42</wall_thermal_conductivity>"
        "</p>");
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> parameters;
    EXPECT_ANY_THROW(BHEns::createPipe(cfg.getConfigSubtree("p"), parameters));
}

/// A zero borehole diameter must be caught by sampleStrictPositive when the
/// assembler asks for cross-section areas or thermal resistances.
TEST(BHEValidation, Bhe1PRejectsZeroDiameter)
{
    ParameterLib::ConstantParameter<double> const zero_diameter{"diameter",
                                                                0.0};
    auto bhe = make1PWith(zero_diameter, kOkWallTc);
    auto const pos = makePos(0);

    EXPECT_THROW(bhe.crossSectionAreas(pos), std::runtime_error);
    EXPECT_THROW(bhe.thermalResistances(pos), std::runtime_error);
}

TEST(BHEValidation, Bhe1PRejectsNegativeDiameter)
{
    ParameterLib::ConstantParameter<double> const neg_diameter{"diameter",
                                                               -0.1};
    auto bhe = make1PWith(neg_diameter, kOkWallTc);
    auto const pos = makePos(0);

    EXPECT_THROW(bhe.crossSectionAreas(pos), std::runtime_error);
}

TEST(BHEValidation, Bhe1PRejectsNaNDiameter)
{
    ParameterLib::ConstantParameter<double> const nan_diameter{
        "diameter", std::numeric_limits<double>::quiet_NaN()};
    auto bhe = make1PWith(nan_diameter, kOkWallTc);
    auto const pos = makePos(0);

    EXPECT_THROW(bhe.crossSectionAreas(pos), std::runtime_error);
}

TEST(BHEValidation, Bhe1PRejectsInfDiameter)
{
    ParameterLib::ConstantParameter<double> const inf_diameter{
        "diameter", std::numeric_limits<double>::infinity()};
    auto bhe = make1PWith(inf_diameter, kOkWallTc);
    auto const pos = makePos(0);

    EXPECT_THROW(bhe.thermalResistances(pos), std::runtime_error);
}

TEST(BHEValidation, Bhe1PRejectsNonPositiveWallTC)
{
    ParameterLib::ConstantParameter<double> const zero_wall_tc{"wall_tc", 0.0};
    auto bhe = make1PWith(kOkDiameter, zero_wall_tc);
    auto const pos = makePos(0);

    EXPECT_THROW(bhe.thermalResistances(pos), std::runtime_error);
}

TEST(BHEValidation, Bhe1PRejectsNaNWallTC)
{
    ParameterLib::ConstantParameter<double> const nan_wall_tc{
        "wall_tc", std::numeric_limits<double>::quiet_NaN()};
    auto bhe = make1PWith(kOkDiameter, nan_wall_tc);
    auto const pos = makePos(0);

    EXPECT_THROW(bhe.thermalResistances(pos), std::runtime_error);
}

/// A borehole diameter smaller than the pipe outside diameter must be rejected
/// by checkBoreholeVsPipeDiameter (1P grout formula requires D > d_outer).
TEST(BHEValidation, Bhe1PRejectsBoreholeSmallerThanPipe)
{
    // Choose D below kPipeOutsideDiameter (0.0458) so the borehole is smaller
    // than the pipe.
    ParameterLib::ConstantParameter<double> const tiny_diameter{"diameter",
                                                                0.04};
    auto bhe = make1PWith(tiny_diameter, kOkWallTc);
    auto const pos = makePos(0);

    EXPECT_THROW(bhe.thermalResistances(pos), std::runtime_error);
}

/// 1U chi formula (Eq. 51) requires D > sqrt(2) * d0; violation must fatal.
TEST(BHEValidation, Bhe1UFailsOnBoreholeVsPipeChiGuard)
{
    // d0 == kPipeOutsideDiameter (0.0458); sqrt(2) * d0 ~ 0.0648, so use
    // D = 0.06 to fail the chi guard while still leaving acosh arguments
    // finite-valued.
    ParameterLib::ConstantParameter<double> const small_diameter{"diameter",
                                                                 0.06};
    auto bhe = make1UWith(small_diameter, kOkWallTc, 0.03);
    auto const pos = makePos(0);

    EXPECT_THROW(bhe.thermalResistances(pos), std::runtime_error);
}

/// 1U R_ar argument (2*s^2 - d0^2) / d0^2 must be > 1, i.e. s > d0.
/// Using s < d0 violates the acosh domain.
TEST(BHEValidation, Bhe1UFailsOnAcoshDomainViolation)
{
    // d0 == kPipeOutsideDiameter (0.0458); choose s = 0.02 so that
    // (2*s*s - d0*d0)/d0*d0 < 1.
    auto bhe = make1UWith(kOkDiameter, kOkWallTc, 0.02);
    auto const pos = makePos(0);

    EXPECT_THROW(bhe.thermalResistances(pos), std::runtime_error);
}

/// Pipes exactly touching (s == d0) make the R_ar acosh argument exactly 1, so
/// acosh(1) = 0 and R_ar = 0. acosh(1) is in-domain but yields a zero
/// resistance and an infinite (1/R) assembly coefficient, so checkAcoshArg
/// must reject arg == 1 strictly.
TEST(BHEValidation, Bhe1UFailsWhenPipesTouch)
{
    // Pipe-center distance == pipe outside diameter => R_ar acosh arg == 1.
    auto bhe = make1UWith(kOkDiameter, kOkWallTc, kPipeOutsideDiameter);
    auto const pos = makePos(0);

    EXPECT_THROW(bhe.thermalResistances(pos), std::runtime_error);
}

/// 2U chi formula (Eq. 38) requires D > 2 * d0.
TEST(BHEValidation, Bhe2UFailsOnBoreholeVsPipeChiGuard)
{
    // 2*d0 ~ 0.0917; use D = 0.08 to violate the guard.
    ParameterLib::ConstantParameter<double> const small_diameter{"diameter",
                                                                 0.08};
    auto bhe = make2UWith(small_diameter, kOkWallTc, 0.03);
    auto const pos = makePos(0);

    EXPECT_THROW(bhe.thermalResistances(pos), std::runtime_error);
}

/// 2U R_ar_1 argument requires s > d0; violate that.
TEST(BHEValidation, Bhe2UFailsOnAcoshDomainViolation)
{
    auto bhe = make2UWith(kOkDiameter, kOkWallTc, 0.02);
    auto const pos = makePos(0);

    EXPECT_THROW(bhe.thermalResistances(pos), std::runtime_error);
}

/// As for 1U: s == d0 makes the R_ar_1 acosh argument exactly 1 and must be
/// rejected strictly.
TEST(BHEValidation, Bhe2UFailsWhenPipesTouch)
{
    auto bhe = make2UWith(kOkDiameter, kOkWallTc, kPipeOutsideDiameter);
    auto const pos = makePos(0);

    EXPECT_THROW(bhe.thermalResistances(pos), std::runtime_error);
}

/// U-type resistance formulas (R_g, R_ar) assume a single pipe outside
/// diameter, so the constructor must reject inlet/outlet pipes whose outside
/// diameters differ. Here the wall thicknesses differ, so do the outside
/// diameters. (Differing wall *conductivity* alone is allowed; see
/// Bhe1UOutletWallConductivityAffectsOutletLeg below.)
TEST(BHEValidation, Bhe1URejectsUnequalPipeOutsideDiameters)
{
    BHEns::BoreholeGeometry const borehole{8.0, kOkDiameter};
    BHEns::Pipe const inlet{kPipeInsideDiameter, kPipeWallThickness, kOkWallTc};
    BHEns::Pipe const outlet{kPipeInsideDiameter, kPipeWallThickness + 0.001,
                             kOkWallTc};
    BHEns::PipeConfigurationUType const pipes{inlet, outlet, 0.06, 0.001};
    auto const control = makeConstantControl(kTemperatureParam, kFlowRateParam);

    EXPECT_THROW(
        (BHEns::BHE_1U{borehole, kRefrigerant, kGrout, control, pipes, false}),
        std::runtime_error);
}

/// As for 1U: a 2U BHE must reject inlet/outlet pipes with different outside
/// diameters.
TEST(BHEValidation, Bhe2URejectsUnequalPipeOutsideDiameters)
{
    BHEns::BoreholeGeometry const borehole{8.0, kOkDiameter};
    BHEns::Pipe const inlet{kPipeInsideDiameter, kPipeWallThickness, kOkWallTc};
    BHEns::Pipe const outlet{kPipeInsideDiameter, kPipeWallThickness + 0.001,
                             kOkWallTc};
    BHEns::PipeConfigurationUType const pipes{inlet, outlet, 0.06, 0.001};
    auto const control = makeConstantControl(kTemperatureParam, kFlowRateParam);

    EXPECT_THROW(
        (BHEns::BHE_2U{borehole, kRefrigerant, kGrout, control, pipes, false}),
        std::runtime_error);
}

/// The inlet leg (R_fig, unknown 0) and outlet leg (R_fog, unknown 1) each use
/// their own pipe wall thermal conductivity. With identical inlet/outlet
/// conductivities the two legs are equal; a distinct outlet conductivity must
/// make R_fog differ from R_fig. Guards against regressing to sampling only the
/// inlet conductivity for both legs.
TEST(BHEThermalResistances, Bhe1UOutletWallConductivityAffectsOutletLeg)
{
    ParameterLib::ConstantParameter<double> const inlet_tc{"inlet_tc", 0.42};
    ParameterLib::ConstantParameter<double> const outlet_tc{"outlet_tc", 4.2};
    auto const pos = makePos(0);

    auto const bhe_equal =
        make1UWithWallTCs(kOkDiameter, inlet_tc, inlet_tc, 0.06);
    auto const R_equal = bhe_equal.thermalResistances(pos);
    EXPECT_DOUBLE_EQ(R_equal[0], R_equal[1]);

    auto const bhe_diff =
        make1UWithWallTCs(kOkDiameter, inlet_tc, outlet_tc, 0.06);
    auto const R_diff = bhe_diff.thermalResistances(pos);
    // Inlet leg unchanged (same inlet conductivity), outlet leg changed.
    EXPECT_DOUBLE_EQ(R_diff[0], R_equal[0]);
    EXPECT_NE(R_diff[1], R_equal[1]);
    EXPECT_LT(R_diff[1], R_diff[0]);  // higher conductivity -> lower resistance
}

TEST(BHEThermalResistances, Bhe2UOutletWallConductivityAffectsOutletLeg)
{
    ParameterLib::ConstantParameter<double> const inlet_tc{"inlet_tc", 0.42};
    ParameterLib::ConstantParameter<double> const outlet_tc{"outlet_tc", 4.2};
    auto const pos = makePos(0);

    auto const bhe_equal =
        make2UWithWallTCs(kOkDiameter, inlet_tc, inlet_tc, 0.06);
    auto const R_equal = bhe_equal.thermalResistances(pos);
    EXPECT_DOUBLE_EQ(R_equal[0], R_equal[1]);

    auto const bhe_diff =
        make2UWithWallTCs(kOkDiameter, inlet_tc, outlet_tc, 0.06);
    auto const R_diff = bhe_diff.thermalResistances(pos);
    EXPECT_DOUBLE_EQ(R_diff[0], R_equal[0]);
    EXPECT_NE(R_diff[1], R_equal[1]);
    EXPECT_LT(R_diff[1], R_diff[0]);
}

namespace
{
/// Build a unit vector tilted `angle_deg` from the canonical vertical-down
/// direction, rotated around +y. angle_deg=0 gives (0, 0, -1), the
/// elem_direction of a vertical mesh with node 0 at the top and node 1 at
/// depth.
Eigen::Vector3d tiltedDownward(double const angle_deg)
{
    double const theta = angle_deg * std::numbers::pi / 180.0;
    return {std::sin(theta), 0.0, -std::cos(theta)};
}
}  // namespace

/// pipeAdvectionVectors must follow the uniform per-leg formula
///     rho_r * Cp_r * flowLegs()[k] * elem_direction
/// for every BHE type. Pre-bhe_param, BHE_1U and BHE_2U ignored
/// elem_direction and hard-coded (0, 0, +/-v), silently making tilted
/// boreholes wrong. This test locks the fix in.

TEST(BHE_1U_Advection, HonorsElemDirection)
{
    auto bhe = make1UWith(kOkDiameter, kOkWallTc, 0.06);

    auto const e_vert = tiltedDownward(0.0);
    auto const e_tilt = tiltedDownward(30.0);

    auto const legs = bhe.flowLegs();
    double const rho_Cp =
        kRefrigerant.density * kRefrigerant.specific_heat_capacity;

    // The two legs of a 1U have signed velocities with opposite signs and
    // equal magnitudes; the active flow must be non-trivial.
    ASSERT_NE(legs[0], 0.0);
    EXPECT_DOUBLE_EQ(legs[0], -legs[1]);

    auto const adv_vert = bhe.pipeAdvectionVectors(e_vert);
    auto const adv_tilt = bhe.pipeAdvectionVectors(e_tilt);

    // Uniform formula: each leg's advection vector equals
    // rho*Cp*v_signed*elem_direction.
    for (int k = 0; k < 2; ++k)
    {
        Eigen::Vector3d const expected_vert = rho_Cp * legs[k] * e_vert;
        Eigen::Vector3d const expected_tilt = rho_Cp * legs[k] * e_tilt;
        EXPECT_NEAR((adv_vert[k] - expected_vert).norm(), 0.0,
                    expected_vert.norm() * 1e-14);
        EXPECT_NEAR((adv_tilt[k] - expected_tilt).norm(), 0.0,
                    expected_tilt.norm() * 1e-14);
    }

    // Grout unknowns (2, 3) carry no advection.
    EXPECT_DOUBLE_EQ(adv_vert[2].norm(), 0.0);
    EXPECT_DOUBLE_EQ(adv_vert[3].norm(), 0.0);
    EXPECT_DOUBLE_EQ(adv_tilt[2].norm(), 0.0);
    EXPECT_DOUBLE_EQ(adv_tilt[3].norm(), 0.0);

    // The pre-fix code hard-coded (0,0,+/-v) and would have returned the
    // same vector for vert and tilt. The fix must yield distinct vectors.
    EXPECT_GT((adv_vert[0] - adv_tilt[0]).norm(), 1e-6);
}

TEST(BHE_2U_Advection, HonorsElemDirection)
{
    auto bhe = make2UWith(kOkDiameter, kOkWallTc, 0.06);

    auto const e_vert = tiltedDownward(0.0);
    auto const e_tilt = tiltedDownward(30.0);

    auto const legs = bhe.flowLegs();
    double const rho_Cp =
        kRefrigerant.density * kRefrigerant.specific_heat_capacity;

    // 2U: { +v, +v, -v, -v } -- two downflow legs and two upflow legs of
    // equal magnitude.
    ASSERT_NE(legs[0], 0.0);
    EXPECT_DOUBLE_EQ(legs[0], legs[1]);
    EXPECT_DOUBLE_EQ(legs[0], -legs[2]);
    EXPECT_DOUBLE_EQ(legs[0], -legs[3]);

    auto const adv_vert = bhe.pipeAdvectionVectors(e_vert);
    auto const adv_tilt = bhe.pipeAdvectionVectors(e_tilt);

    for (int k = 0; k < 4; ++k)
    {
        Eigen::Vector3d const expected_vert = rho_Cp * legs[k] * e_vert;
        Eigen::Vector3d const expected_tilt = rho_Cp * legs[k] * e_tilt;
        EXPECT_NEAR((adv_vert[k] - expected_vert).norm(), 0.0,
                    expected_vert.norm() * 1e-14);
        EXPECT_NEAR((adv_tilt[k] - expected_tilt).norm(), 0.0,
                    expected_tilt.norm() * 1e-14);
    }

    // Grout unknowns (4..7) carry no advection.
    for (int g = 4; g < 8; ++g)
    {
        EXPECT_DOUBLE_EQ(adv_vert[g].norm(), 0.0);
        EXPECT_DOUBLE_EQ(adv_tilt[g].norm(), 0.0);
    }

    // Tilted result must differ from the vertical one.
    EXPECT_GT((adv_vert[0] - adv_tilt[0]).norm(), 1e-6);
}

TEST(BHE_1P_Advection, HonorsElemDirection)
{
    auto bhe = make1PWith(kOkDiameter, kOkWallTc);

    auto const e_vert = tiltedDownward(0.0);
    auto const e_tilt = tiltedDownward(30.0);

    auto const legs = bhe.flowLegs();
    double const rho_Cp =
        kRefrigerant.density * kRefrigerant.specific_heat_capacity;

    ASSERT_NE(legs[0], 0.0);

    auto const adv_vert = bhe.pipeAdvectionVectors(e_vert);
    auto const adv_tilt = bhe.pipeAdvectionVectors(e_tilt);

    Eigen::Vector3d const expected_vert = rho_Cp * legs[0] * e_vert;
    Eigen::Vector3d const expected_tilt = rho_Cp * legs[0] * e_tilt;
    EXPECT_NEAR((adv_vert[0] - expected_vert).norm(), 0.0,
                expected_vert.norm() * 1e-14);
    EXPECT_NEAR((adv_tilt[0] - expected_tilt).norm(), 0.0,
                expected_tilt.norm() * 1e-14);

    // Grout unknown (1) carries no advection.
    EXPECT_DOUBLE_EQ(adv_vert[1].norm(), 0.0);
    EXPECT_DOUBLE_EQ(adv_tilt[1].norm(), 0.0);
}

namespace
{
/// Standard coaxial pipe configuration shared by the CXA/CXC coverage tests.
BHEns::PipeConfigurationCoaxial makeCoaxialPipes()
{
    static ParameterLib::ConstantParameter<double> const inner_tc{"inner_tc",
                                                                  0.001};
    static ParameterLib::ConstantParameter<double> const outer_tc{"outer_tc",
                                                                  1.3};
    BHEns::Pipe const inner_pipe{0.09532, 0.00734, inner_tc};
    BHEns::Pipe const outer_pipe{0.16626, 0.00587, outer_tc};
    return BHEns::PipeConfigurationCoaxial{inner_pipe, outer_pipe, 0.001};
}

/// Exercise the per-element property accessors that the assembler calls for
/// coaxial BHEs (thermalResistances, pipeHeatConductions, pipeAdvectionVectors)
/// so both CXA and CXC have unit coverage for them.
template <typename CoaxialBHE>
void checkCoaxialElementWiseProperties()
{
    BHEns::BoreholeGeometry const borehole{8.0, kOkDiameter};
    auto const control = makeConstantControl(kTemperatureParam, kFlowRateParam);
    CoaxialBHE const bhe{borehole, kRefrigerant,       kGrout,
                         control,  makeCoaxialPipes(), false};
    auto const pos = makePos(0);

    // thermalResistances: 3 finite, positive values.
    auto const R = bhe.thermalResistances(pos);
    ASSERT_EQ(R.size(), 3u);
    for (int u = 0; u < 3; ++u)
    {
        EXPECT_TRUE(std::isfinite(R[u])) << "unknown=" << u;
        EXPECT_GT(R[u], 0.0) << "unknown=" << u;
    }

    // pipeHeatConductions: 3 finite values.
    auto const lambda = bhe.pipeHeatConductions();
    for (int u = 0; u < 3; ++u)
    {
        EXPECT_TRUE(std::isfinite(lambda[u])) << "unknown=" << u;
    }

    // pipeAdvectionVectors: each fluid leg follows the uniform per-leg formula
    // rho_r * Cp_r * flowLegs()[k] * elem_direction; the grout unknown (2)
    // carries no advection.
    auto const e_tilt = tiltedDownward(30.0);
    auto const legs = bhe.flowLegs();
    double const rho_Cp =
        kRefrigerant.density * kRefrigerant.specific_heat_capacity;
    auto const adv = bhe.pipeAdvectionVectors(e_tilt);
    for (int k = 0; k < 2; ++k)
    {
        Eigen::Vector3d const expected = rho_Cp * legs[k] * e_tilt;
        EXPECT_NEAR((adv[k] - expected).norm(), 0.0, expected.norm() * 1e-14)
            << "leg=" << k;
    }
    EXPECT_DOUBLE_EQ(adv[2].norm(), 0.0);
}
}  // namespace

TEST(BHECommonCoaxial, CxaElementWiseProperties)
{
    checkCoaxialElementWiseProperties<BHEns::BHE_CXA>();
}

TEST(BHECommonCoaxial, CxcElementWiseProperties)
{
    checkCoaxialElementWiseProperties<BHEns::BHE_CXC>();
}

namespace
{
/// A self-contained line-element chain: owns its nodes and elements so the
/// test can hand the raw `Element*` pointers to
/// `findBHEEndpointsFromElementOrdering` without ever building a `Mesh`.
struct LineChain
{
    std::vector<std::unique_ptr<MeshLib::Node>> nodes;
    std::vector<std::unique_ptr<MeshLib::Element>> elements;
    std::vector<MeshLib::Element*> raw_elements;
};

/// Build a vertical line chain of `n_elements` elements. Node `k` sits at
/// `z = -k`, and element `k` runs from node `k` (local 0) to node `k+1`
/// (local 1) -- i.e. the canonical "node 0 at higher z, node 1 at lower z"
/// convention used by every existing BHE benchmark mesh.
LineChain makeVerticalLineChain(std::size_t n_elements)
{
    LineChain c;
    for (std::size_t i = 0; i <= n_elements; ++i)
    {
        c.nodes.push_back(std::make_unique<MeshLib::Node>(
            std::array<double, 3>{0.0, 0.0, -static_cast<double>(i)}, i));
    }
    for (std::size_t ei = 0; ei < n_elements; ++ei)
    {
        std::array<MeshLib::Node*, 2> en{c.nodes[ei].get(),
                                         c.nodes[ei + 1].get()};
        c.elements.push_back(std::make_unique<MeshLib::Line>(en, ei));
        c.raw_elements.push_back(c.elements.back().get());
    }
    return c;
}
}  // namespace

TEST(BHEEndpoints, SingleElementChainPicksNode0AsInlet)
{
    auto c = makeVerticalLineChain(1);
    auto const ep =
        ProcessLib::HeatTransportBHE::findBHEEndpointsFromElementOrdering(
            c.raw_elements);
    EXPECT_EQ(ep.inlet, c.nodes[0].get());
    EXPECT_EQ(ep.outlet, c.nodes[1].get());
}

TEST(BHEEndpoints, MultiElementChainIdentifiesEndpointsRegardlessOfVectorOrder)
{
    auto c = makeVerticalLineChain(5);
    // The helper must walk the chain via node adjacency, not by trusting
    // the order in which the elements appear in the input vector.
    std::reverse(c.raw_elements.begin(), c.raw_elements.end());
    auto const ep =
        ProcessLib::HeatTransportBHE::findBHEEndpointsFromElementOrdering(
            c.raw_elements);
    EXPECT_EQ(ep.inlet, c.nodes[0].get());
    EXPECT_EQ(ep.outlet, c.nodes[5].get());
}

TEST(BHEEndpoints, RejectsReversedInteriorElement)
{
    auto c = makeVerticalLineChain(3);
    // Replace the middle element (originally nodes[1] -> nodes[2]) with a
    // reversed copy (nodes[2] -> nodes[1]). The chain walk should refuse
    // to advance from nodes[1] because no unused element has nodes[1] as
    // its local node 0.
    std::array<MeshLib::Node*, 2> reversed{c.nodes[2].get(), c.nodes[1].get()};
    c.elements[1] = std::make_unique<MeshLib::Line>(reversed, 1);
    c.raw_elements[1] = c.elements[1].get();
    EXPECT_THROW(
        ProcessLib::HeatTransportBHE::findBHEEndpointsFromElementOrdering(
            c.raw_elements),
        std::runtime_error);
}

TEST(BHEEndpoints, RejectsBranchedChain)
{
    auto c = makeVerticalLineChain(2);  // nodes 0,1,2 with elements 0-1, 1-2
    // Attach a third element that branches off from node 1 sideways.
    c.nodes.push_back(std::make_unique<MeshLib::Node>(
        std::array<double, 3>{1.0, 0.0, -1.0}, 3));
    std::array<MeshLib::Node*, 2> branch{c.nodes[1].get(), c.nodes[3].get()};
    c.elements.push_back(std::make_unique<MeshLib::Line>(branch, 2));
    c.raw_elements.push_back(c.elements.back().get());
    EXPECT_THROW(
        ProcessLib::HeatTransportBHE::findBHEEndpointsFromElementOrdering(
            c.raw_elements),
        std::runtime_error);
}

TEST(BHEEndpoints, RejectsDisconnectedChain)
{
    // Two separate single-element chains with no shared node -- four
    // endpoints in total, not two, so the helper must reject.
    LineChain c;
    for (std::size_t i = 0; i < 4; ++i)
    {
        c.nodes.push_back(std::make_unique<MeshLib::Node>(
            std::array<double, 3>{0.0, 0.0, -static_cast<double>(i)}, i));
    }
    std::array<MeshLib::Node*, 2> e0{c.nodes[0].get(), c.nodes[1].get()};
    std::array<MeshLib::Node*, 2> e1{c.nodes[2].get(), c.nodes[3].get()};
    c.elements.push_back(std::make_unique<MeshLib::Line>(e0, 0));
    c.elements.push_back(std::make_unique<MeshLib::Line>(e1, 1));
    c.raw_elements.push_back(c.elements[0].get());
    c.raw_elements.push_back(c.elements[1].get());
    EXPECT_THROW(
        ProcessLib::HeatTransportBHE::findBHEEndpointsFromElementOrdering(
            c.raw_elements),
        std::runtime_error);
}

TEST(BHEEndpoints, HorizontalBheUsesElementOrderingForEndpoints)
{
    // A purely horizontal chain: every node shares z = 0, so the retired
    // z-coordinate wellhead heuristic could not decide which end is the
    // inlet. Endpoint selection now follows the line-element node ordering,
    // so a horizontal BHE is unambiguous as long as its elements are
    // consistently oriented (local node 0 -> local node 1 along the chain).
    LineChain c;
    std::size_t const n_elements = 4;
    for (std::size_t i = 0; i <= n_elements; ++i)
    {
        c.nodes.push_back(std::make_unique<MeshLib::Node>(
            std::array<double, 3>{static_cast<double>(i), 0.0, 0.0}, i));
    }
    for (std::size_t ei = 0; ei < n_elements; ++ei)
    {
        std::array<MeshLib::Node*, 2> en{c.nodes[ei].get(),
                                         c.nodes[ei + 1].get()};
        c.elements.push_back(std::make_unique<MeshLib::Line>(en, ei));
        c.raw_elements.push_back(c.elements.back().get());
    }
    auto const ep =
        ProcessLib::HeatTransportBHE::findBHEEndpointsFromElementOrdering(
            c.raw_elements);
    EXPECT_EQ(ep.inlet, c.nodes[0].get());
    EXPECT_EQ(ep.outlet, c.nodes[n_elements].get());
}

/// getBHEDataInMesh must accept a horizontal BHE (all nodes share z) and
/// extract its line elements and nodes as one material group. The removed
/// ProcessLibBHEMeshUtils.HorizontalBheWithAmbiguousWellheadFails asserted the
/// opposite (a throw), which became invalid once endpoint selection moved from
/// a z-coordinate wellhead to element node ordering
/// (findBHEEndpointsFromElementOrdering); horizontal geometry is now valid.
TEST(ProcessLibBHEMeshUtils, HorizontalBheIsExtractedFromMesh)
{
    // MeshLib::Mesh takes ownership of the raw Node* / Element* pointers.
    std::vector<MeshLib::Node*> nodes{new MeshLib::Node(0.0, 0.0, 0.0, 0),
                                      new MeshLib::Node(1.0, 0.0, 0.0, 1),
                                      new MeshLib::Node(2.0, 0.0, 0.0, 2)};
    std::vector<MeshLib::Element*> elements{
        new MeshLib::Line(std::array<MeshLib::Node*, 2>{{nodes[0], nodes[1]}}),
        new MeshLib::Line(std::array<MeshLib::Node*, 2>{{nodes[1], nodes[2]}})};

    auto mesh = std::make_unique<MeshLib::Mesh>(
        "horizontal_bhe", nodes, elements, true /*compute_element_neighbors*/);
    MeshLib::addPropertyToMesh<int>(*mesh, "MaterialIDs",
                                    MeshLib::MeshItemType::Cell, 1, {{0, 0}});

    // Was expected to throw (ambiguous wellhead); horizontal is now valid.
    auto const data = ProcessLib::HeatTransportBHE::getBHEDataInMesh(*mesh);

    ASSERT_EQ(data.BHE_mat_IDs.size(), 1u);
    EXPECT_EQ(data.BHE_mat_IDs.front(), 0);
    ASSERT_EQ(data.BHE_elements.size(), 1u);
    EXPECT_EQ(data.BHE_elements.front().size(), 2u);
    ASSERT_EQ(data.BHE_nodes.size(), 1u);
    EXPECT_EQ(data.BHE_nodes.front().size(), 3u);
}
