// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <memory>
#include <vector>

#include "BaseLib/ConfigTree.h"
#include "BaseLib/Error.h"
#include "MeshLib/Elements/Line.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"
#include "MeshLib/Utils/addPropertyToMesh.h"
#include "ParameterLib/ConstantParameter.h"
#include "ParameterLib/FunctionParameter.h"
#include "ParameterLib/Parameter.h"
#include "ParameterLib/SpatialPosition.h"
#include "ProcessLib/HeatTransportBHE/BHE/BHESectionUtils.h"
#include "ProcessLib/HeatTransportBHE/BHE/BHE_1P.h"
#include "ProcessLib/HeatTransportBHE/BHE/BHE_1U.h"
#include "ProcessLib/HeatTransportBHE/BHE/BHE_CXA.h"
#include "ProcessLib/HeatTransportBHE/BHE/BoreholeGeometry.h"
#include "ProcessLib/HeatTransportBHE/BHE/FlowAndTemperatureControl.h"
#include "ProcessLib/HeatTransportBHE/BHE/GroutParameters.h"
#include "ProcessLib/HeatTransportBHE/BHE/MeshUtils.h"
#include "ProcessLib/HeatTransportBHE/BHE/Pipe.h"
#include "ProcessLib/HeatTransportBHE/BHE/PipeConfigurationCoaxial.h"
#include "ProcessLib/HeatTransportBHE/BHE/PipeConfigurationUType.h"
#include "ProcessLib/HeatTransportBHE/BHE/RefrigerantProperties.h"
#include "Tests/TestTools.h"

namespace
{
BaseLib::ConfigTree makeConfigTree(char const* xml)
{
    auto ptree = Tests::readXml(xml);
    return BaseLib::ConfigTree(std::move(ptree), "TestBHEElementWiseProperties",
                               BaseLib::ConfigTree::onerror,
                               BaseLib::ConfigTree::onwarning);
}

std::vector<std::unique_ptr<MeshLib::Node>> createInclinedBheNodes()
{
    std::vector<std::unique_ptr<MeshLib::Node>> nodes;
    nodes.emplace_back(std::make_unique<MeshLib::Node>(0.0, 0.0, 0.0, 0));
    nodes.emplace_back(std::make_unique<MeshLib::Node>(1.0, 0.0, -1.0, 1));
    nodes.emplace_back(std::make_unique<MeshLib::Node>(2.0, 0.0, -2.0, 2));
    return nodes;
}

std::vector<std::unique_ptr<MeshLib::Node>> createVerticalBheNodesAtX(
    double const x, std::size_t const id_offset)
{
    std::vector<std::unique_ptr<MeshLib::Node>> nodes;
    nodes.emplace_back(std::make_unique<MeshLib::Node>(x, 0.0, 0.0, id_offset));
    nodes.emplace_back(
        std::make_unique<MeshLib::Node>(x, 0.0, -4.0, id_offset + 1));
    nodes.emplace_back(
        std::make_unique<MeshLib::Node>(x, 0.0, -8.0, id_offset + 2));
    return nodes;
}

std::vector<MeshLib::Node*> toNodePtrs(
    std::vector<std::unique_ptr<MeshLib::Node>> const& nodes)
{
    std::vector<MeshLib::Node*> ptrs;
    ptrs.reserve(nodes.size());
    for (auto const& node : nodes)
    {
        ptrs.push_back(node.get());
    }
    return ptrs;
}

/// Sort BHE nodes by z-coordinate descending (wellhead = top = largest z).
/// Test-local helper; not needed in production code.
std::vector<MeshLib::Node*> sortNodesByDepth(
    std::vector<MeshLib::Node*> const& bhe_nodes)
{
    if (bhe_nodes.empty())
    {
        OGS_FATAL(
            "BHE node list is empty. Cannot evaluate diameter parameter.");
    }

    std::vector<MeshLib::Node*> sorted_nodes(bhe_nodes.begin(),
                                             bhe_nodes.end());
    std::sort(sorted_nodes.begin(), sorted_nodes.end(),
              [](MeshLib::Node const* const a, MeshLib::Node const* const b)
              { return (*a)[2] > (*b)[2]; });
    return sorted_nodes;
}

auto makeProfile(ParameterLib::Parameter<double> const& parameter,
                 std::vector<MeshLib::Node*> const& nodes)
{
    using namespace ProcessLib::HeatTransportBHE::BHE;
    auto const sorted = sortNodesByDepth(nodes);
    auto const distances = cumulativeDistances(sorted);
    std::vector<double> diameters;
    diameters.reserve(sorted.size());
    for (auto const* node : sorted)
    {
        ParameterLib::SpatialPosition pos;
        pos.setNodeID(node->getID());
        pos.setCoordinates(*node);
        diameters.push_back(parameter(0.0, pos)[0]);
    }
    return groupSections(distances, diameters);
}
}  // namespace

TEST(ProcessLibBHEMeshUtils, ArcLengthDistanceForInclinedBhe)
{
    // MeshLib::Mesh takes ownership of raw Node* and Element* pointers.
    std::vector<MeshLib::Node*> nodes;
    nodes.push_back(new MeshLib::Node(0.0, 0.0, 0.0, 0));
    nodes.push_back(new MeshLib::Node(1.0, 0.0, -1.0, 1));
    nodes.push_back(new MeshLib::Node(2.0, 0.0, -2.0, 2));

    std::vector<MeshLib::Element*> elements;
    elements.push_back(
        new MeshLib::Line(std::array<MeshLib::Node*, 2>{{nodes[0], nodes[1]}}));
    elements.push_back(
        new MeshLib::Line(std::array<MeshLib::Node*, 2>{{nodes[1], nodes[2]}}));

    auto mesh = std::make_unique<MeshLib::Mesh>(
        "inclined_bhe", nodes, elements, true /*compute_element_neighbors*/);

    MeshLib::addPropertyToMesh<int>(*mesh, "MaterialIDs",
                                    MeshLib::MeshItemType::Cell, 1, {{0, 0}});

    auto const bhe_mesh_data =
        ProcessLib::HeatTransportBHE::getBHEDataInMesh(*mesh);

    double const segment_length = std::sqrt(2.0);
    EXPECT_NEAR(bhe_mesh_data.BHE_element_distances_from_wellhead.at(
                    elements[0]->getID()),
                0.5 * segment_length, 1e-12);
    EXPECT_NEAR(bhe_mesh_data.BHE_element_distances_from_wellhead.at(
                    elements[1]->getID()),
                1.5 * segment_length, 1e-12);
}

TEST(ProcessLibBHEMeshUtils, HorizontalBheWithAmbiguousWellheadFails)
{
    // Use node IDs such that the interior node has the smallest id.
    // Both endpoints share the same z-coordinate, so the wellhead choice is
    // ambiguous and must fail.
    // MeshLib::Mesh takes ownership of raw Node* and Element* pointers.
    std::vector<MeshLib::Node*> nodes;
    nodes.push_back(new MeshLib::Node(0.0, 0.0, 0.0, 1));  // endpoint
    nodes.push_back(new MeshLib::Node(1.0, 0.0, 0.0, 0));  // interior
    nodes.push_back(new MeshLib::Node(2.0, 0.0, 0.0, 2));  // endpoint

    std::vector<MeshLib::Element*> elements;
    elements.push_back(
        new MeshLib::Line(std::array<MeshLib::Node*, 2>{{nodes[0], nodes[1]}}));
    elements.push_back(
        new MeshLib::Line(std::array<MeshLib::Node*, 2>{{nodes[1], nodes[2]}}));

    auto mesh = std::make_unique<MeshLib::Mesh>(
        "horizontal_bhe", nodes, elements, true /*compute_element_neighbors*/);

    MeshLib::addPropertyToMesh<int>(*mesh, "MaterialIDs",
                                    MeshLib::MeshItemType::Cell, 1, {{0, 0}});

    EXPECT_ANY_THROW(ProcessLib::HeatTransportBHE::getBHEDataInMesh(*mesh));
}

TEST(ProcessLibBHEGeometry, InlineDiameterCreatesParameterAndSingleSection)
{
    auto config = makeConfigTree(
        "<borehole>"
        "  <length>20.0</length>"
        "  <diameter>0.16</diameter>"
        "</borehole>");

    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> parameters;
    auto owned_nodes = createInclinedBheNodes();
    auto const bhe_nodes = toNodePtrs(owned_nodes);

    auto const borehole_geometry =
        ProcessLib::HeatTransportBHE::BHE::createBoreholeGeometry(
            config.getConfigSubtree("borehole"), parameters, bhe_nodes);

    EXPECT_EQ(1, borehole_geometry.sections.getNumberOfSections());
    EXPECT_DOUBLE_EQ(0.16, borehole_geometry.sections.section_diameters[0]);
    EXPECT_EQ(1u, parameters.size());
}

TEST(ProcessLibBHEPipe, CreatePipeRejectsInvalidParameters)
{
    auto config_bad_diameter = makeConfigTree(
        "<p>"
        "  <diameter>-0.04</diameter>"
        "  <wall_thickness>0.0029</wall_thickness>"
        "  <wall_thermal_conductivity>0.42</wall_thermal_conductivity>"
        "</p>");
    EXPECT_ANY_THROW(ProcessLib::HeatTransportBHE::BHE::createPipe(
        config_bad_diameter.getConfigSubtree("p")));

    auto config_bad_thickness = makeConfigTree(
        "<p>"
        "  <diameter>0.04</diameter>"
        "  <wall_thickness>-0.001</wall_thickness>"
        "  <wall_thermal_conductivity>0.42</wall_thermal_conductivity>"
        "</p>");
    EXPECT_ANY_THROW(ProcessLib::HeatTransportBHE::BHE::createPipe(
        config_bad_thickness.getConfigSubtree("p")));

    auto config_bad_conductivity = makeConfigTree(
        "<p>"
        "  <diameter>0.04</diameter>"
        "  <wall_thickness>0.0029</wall_thickness>"
        "  <wall_thermal_conductivity>0.0</wall_thermal_conductivity>"
        "</p>");
    EXPECT_ANY_THROW(ProcessLib::HeatTransportBHE::BHE::createPipe(
        config_bad_conductivity.getConfigSubtree("p")));
}

/// Regression test: CXA BHE with a sectioned borehole (2 sections) and
/// constant-diameter pipes must produce finite thermal resistances for all
/// borehole sections.
TEST(BHECommonCoaxial, CxaBoreholeOnlySectionsDoesNotCrash)
{
    using namespace ProcessLib::HeatTransportBHE::BHE;

    // Borehole: 2 sections (mirrors "if(z > -4, 0.22, 0.20)" function param)
    BoreholeGeometry const borehole{8.0, {{0.0, 4.0}, {0.22, 0.20}}};

    RefrigerantProperties const refrigerant{0.00067418, 992.92, 0.62863, 4198.0,
                                            25.0};

    GroutParameters const grout{2190.0, 0.0, 1735.16, 0.806};

    // Pipes: constant scalar diameters
    Pipe const inner_pipe{0.09532, 0.00734, 0.001};
    Pipe const outer_pipe{0.16626, 0.00587, 1.3};
    PipeConfigurationCoaxial const pipes{inner_pipe, outer_pipe, 0.001};

    // InflowTemperature control requires Parameter<double> references.
    // ParameterLib::ConstantParameter supplies fixed scalar values.
    ParameterLib::ConstantParameter<double> const flow_rate_param{"flow_rate",
                                                                  2.0e-4};
    ParameterLib::ConstantParameter<double> const temperature_param{
        "temperature", 25.0};
    FlowAndTemperatureControl const control =
        InflowTemperature{temperature_param, flow_rate_param, 0.0};

    // Construction calls updateHeatTransferCoefficients, which calls
    // calcThermalResistances for each borehole section. Without the fix this
    // crashes; with the fix it completes normally.
    BHE_CXA const bhe{borehole, refrigerant, grout, control, pipes, false};

    EXPECT_EQ(2, bhe.getNumberOfSections());
    // Both sections must yield finite thermal resistance values.
    EXPECT_TRUE(std::isfinite(bhe.thermalResistanceAtSection(0, 0)));
    EXPECT_TRUE(std::isfinite(bhe.thermalResistanceAtSection(0, 1)));

    // Grout cross-section area (index 2) must differ between sections
    // because the borehole diameters differ (0.22 vs 0.20).
    auto const areas_s0 = bhe.crossSectionAreas(0);
    auto const areas_s1 = bhe.crossSectionAreas(1);
    EXPECT_NE(areas_s0[2], areas_s1[2]);
    // Both grout areas must be positive.
    EXPECT_GT(areas_s0[2], 0.0);
    EXPECT_GT(areas_s1[2], 0.0);
}

/// Regression test: BHE_1U with a sectioned borehole (2 sections) must
/// produce finite thermal resistances for all sections. This exercises
/// the U-type thermal resistance path (thermalResistancesGroutSoil with
/// the chi correction loop).
TEST(BHECommonUType, Bhe1UBoreholeOnlySectionsDoesNotCrash)
{
    using namespace ProcessLib::HeatTransportBHE::BHE;

    // Borehole: 2 sections with different diameters
    BoreholeGeometry const borehole{8.0, {{0.0, 4.0}, {0.22, 0.20}}};

    RefrigerantProperties const refrigerant{0.00067418, 992.92, 0.62863, 4198.0,
                                            25.0};

    GroutParameters const grout{2190.0, 0.0, 1735.16, 0.806};

    // U-type pipes: inlet and outlet with distance between them
    Pipe const inlet{0.04, 0.0029, 0.42};
    Pipe const outlet{0.04, 0.0029, 0.42};
    PipeConfigurationUType const pipes{inlet, outlet, 0.06, 0.001};

    ParameterLib::ConstantParameter<double> const flow_rate_param{"flow_rate",
                                                                  2.0e-4};
    ParameterLib::ConstantParameter<double> const temperature_param{
        "temperature", 25.0};
    FlowAndTemperatureControl const control =
        InflowTemperature{temperature_param, flow_rate_param, 0.0};

    BHE_1U const bhe{borehole, refrigerant, grout, control, pipes, false};

    EXPECT_EQ(2, bhe.getNumberOfSections());
    // All 4 resistance unknowns must be finite for both sections.
    for (int s = 0; s < 2; ++s)
    {
        for (int u = 0; u < 4; ++u)
        {
            EXPECT_TRUE(std::isfinite(bhe.thermalResistanceAtSection(u, s)))
                << "section=" << s << " unknown=" << u;
        }
    }

    // Grout cross-section areas (indices 2,3) must differ between sections.
    auto const areas_s0 = bhe.crossSectionAreas(0);
    auto const areas_s1 = bhe.crossSectionAreas(1);
    EXPECT_NE(areas_s0[2], areas_s1[2]);
    EXPECT_GT(areas_s0[2], 0.0);
    EXPECT_GT(areas_s1[2], 0.0);
}

/// Regression test for grouped BHE definitions (id="*" or multi-id):
/// When the diameter parameter is spatially varying, each BHE must get its own
/// DiameterProfile based on its own node positions.  The grouped BHE code path
/// must not reuse the first BHE's profile for all subsequent IDs.
TEST(ProcessLibBHEGeometry,
     SpatiallyVaryingDiameterProducesDifferentProfilesPerNodeSet)
{
    auto owned_a = createVerticalBheNodesAtX(0.0, 0);
    auto owned_b = createVerticalBheNodesAtX(1.0, 3);
    auto const nodes_a = toNodePtrs(owned_a);
    auto const nodes_b = toNodePtrs(owned_b);

    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>>
        curves;
    ParameterLib::FunctionParameter<double> param(
        "x_dependent_diameter", {"if(x >= 0.5, 0.14, 0.11)"}, curves);

    auto const [bounds_a, diams_a] = makeProfile(param, nodes_a);
    auto const [bounds_b, diams_b] = makeProfile(param, nodes_b);

    // Both should be single-section profiles (all nodes at constant x).
    ASSERT_EQ(1u, bounds_a.size());
    ASSERT_EQ(1u, bounds_b.size());
    ASSERT_EQ(1u, diams_a.size());
    ASSERT_EQ(1u, diams_b.size());
    EXPECT_DOUBLE_EQ(0.0, bounds_a[0]);
    EXPECT_DOUBLE_EQ(0.0, bounds_b[0]);

    // BHE 0 at x=0 should get small diameter.
    EXPECT_DOUBLE_EQ(0.11, diams_a[0]);
    // BHE 1 at x=1 should get large diameter.
    EXPECT_DOUBLE_EQ(0.14, diams_b[0]);

    // Key assertion: profiles must differ.  The grouped-BHE bug incorrectly
    // shares the first ID's profile with all subsequent IDs.
    EXPECT_NE(diams_a[0], diams_b[0]);
}

TEST(ProcessLibBHEGeometry, GroupedBheRebuildForNodesProducesDifferentGeometry)
{
    using namespace ProcessLib::HeatTransportBHE::BHE;

    auto owned_a = createVerticalBheNodesAtX(0.0, 0);
    auto owned_b = createVerticalBheNodesAtX(1.0, 3);
    auto const nodes_a = toNodePtrs(owned_a);
    auto const nodes_b = toNodePtrs(owned_b);

    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>>
        curves;
    ParameterLib::FunctionParameter<double> param(
        "x_dependent_diameter", {"if(x >= 0.5, 0.14, 0.11)"}, curves);

    auto const [bounds_a, diams_a] = makeProfile(param, nodes_a);
    BoreholeGeometry const geometry_a{
        8.0,
        {bounds_a, diams_a},
        static_cast<ParameterLib::Parameter<double> const*>(&param)};

    auto const geometry_b = geometry_a.rebuildForNodes(nodes_b);

    ASSERT_EQ(1, geometry_a.sections.getNumberOfSections());
    ASSERT_EQ(1, geometry_b.sections.getNumberOfSections());
    EXPECT_DOUBLE_EQ(0.0, geometry_a.sections.section_boundaries[0]);
    EXPECT_DOUBLE_EQ(0.0, geometry_b.sections.section_boundaries[0]);
    EXPECT_DOUBLE_EQ(0.11, geometry_a.sections.section_diameters[0]);
    EXPECT_DOUBLE_EQ(0.14, geometry_b.sections.section_diameters[0]);
    EXPECT_NE(geometry_a.sections.section_diameters[0],
              geometry_b.sections.section_diameters[0]);
}

/// Regression test: BHE_1P with a sectioned borehole (2 sections) must
/// produce finite thermal resistances for all sections.
TEST(BHECommon1P, Bhe1PBoreholeOnlySectionsDoesNotCrash)
{
    using namespace ProcessLib::HeatTransportBHE::BHE;

    // Borehole: 2 sections with different diameters
    BoreholeGeometry const borehole{8.0, {{0.0, 4.0}, {0.22, 0.20}}};

    RefrigerantProperties const refrigerant{0.00067418, 992.92, 0.62863, 4198.0,
                                            25.0};

    GroutParameters const grout{2190.0, 0.0, 1735.16, 0.806};

    // Single pipe
    Pipe const single_pipe{0.04, 0.0029, 0.42};
    PipeConfiguration1PType const pipes{single_pipe, 0.001};

    ParameterLib::ConstantParameter<double> const flow_rate_param{"flow_rate",
                                                                  2.0e-4};
    ParameterLib::ConstantParameter<double> const temperature_param{
        "temperature", 25.0};
    FlowAndTemperatureControl const control =
        InflowTemperature{temperature_param, flow_rate_param, 0.0};

    BHE_1P const bhe{borehole, refrigerant, grout, control, pipes, false};

    EXPECT_EQ(2, bhe.getNumberOfSections());
    // Both sections must yield finite thermal resistance values.
    for (int s = 0; s < 2; ++s)
    {
        for (int u = 0; u < 2; ++u)
        {
            EXPECT_TRUE(std::isfinite(bhe.thermalResistanceAtSection(u, s)))
                << "section=" << s << " unknown=" << u;
        }
    }

    // Grout cross-section area (index 1) must differ between sections
    // because the borehole diameters differ (0.22 vs 0.20).
    auto const areas_s0 = bhe.crossSectionAreas(0);
    auto const areas_s1 = bhe.crossSectionAreas(1);
    EXPECT_NE(areas_s0[1], areas_s1[1]);
    EXPECT_GT(areas_s0[1], 0.0);
    EXPECT_GT(areas_s1[1], 0.0);
}
