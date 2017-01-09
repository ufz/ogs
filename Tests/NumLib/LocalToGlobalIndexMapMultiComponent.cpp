/*
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <gtest/gtest.h>

#include <memory>

#include "NumLib/DOF/LocalToGlobalIndexMap.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/Location.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/MeshSubsets.h"
#include "MeshLib/MeshGenerators/MeshGenerator.h"
#include "MeshLib/MeshSearch/NodeSearch.h"
#include "GeoLib/Polyline.h"
#include "GeoLib/GEOObjects.h"

#include "MeshGeoToolsLib/MeshNodeSearcher.h"
#include "MeshGeoToolsLib/BoundaryElementsSearcher.h"

namespace NL = NumLib;
namespace MeL = MeshLib;
namespace MGTL = MeshGeoToolsLib;

class NumLibLocalToGlobalIndexMapMultiDOFTest : public ::testing::Test
{
public:
    static const std::size_t mesh_subdivs = 4;
    NumLibLocalToGlobalIndexMapMultiDOFTest()
    {
        mesh.reset(MeL::MeshGenerator::generateRegularQuadMesh(2.0, mesh_subdivs));
        mesh_items_all_nodes.reset(new MeL::MeshSubset(*mesh, &mesh->getNodes()));

        std::unique_ptr<std::vector<GeoLib::Point*>> ply_pnts(
                    new std::vector<GeoLib::Point*>);
        ply_pnts->push_back(new GeoLib::Point(0.0, 0.0, 0.0));
        ply_pnts->push_back(new GeoLib::Point(1.0, 0.0, 0.0));

        std::string geometry_0("GeometryWithPntsAndPolyline");
        geo_objs.addPointVec(std::move(ply_pnts), geometry_0, nullptr);

        auto ply = new GeoLib::Polyline(*geo_objs.getPointVec(geometry_0));

        ply->addPoint(0);
        ply->addPoint(1);

        std::unique_ptr<std::vector<GeoLib::Polyline*>> plys(
                    new std::vector<GeoLib::Polyline*>);
        plys->push_back(ply);

        geo_objs.addPolylineVec(std::move(plys), geometry_0, nullptr);

        MGTL::MeshNodeSearcher& searcher_nodes = MGTL::MeshNodeSearcher::getMeshNodeSearcher(*mesh);
        MGTL::BoundaryElementsSearcher searcher_elements(*mesh, searcher_nodes);

        auto elems = searcher_elements.getBoundaryElements(*ply);

        // deep copy because the searcher destroys the elements.
        std::transform(elems.cbegin(), elems.cend(),
                       std::back_inserter(boundary_elements),
                       [](MeL::Element const* e) { return e->clone(); });

        std::vector<MeL::Node*> nodes = MeL::getUniqueNodes(boundary_elements);

        mesh_items_boundary.reset(
            mesh_items_all_nodes->getIntersectionByNodes(nodes));
    }

    ~NumLibLocalToGlobalIndexMapMultiDOFTest()
    {
        for (auto e : boundary_elements)
            delete e;
    }

    void initComponents(const unsigned num_components, const unsigned selected_component,
                        const NL::ComponentOrder order)
    {
        assert(selected_component < num_components);

        std::vector<std::unique_ptr<MeshLib::MeshSubsets>> components;
        for (unsigned i=0; i<num_components; ++i)
        {
            components.emplace_back(
                new MeL::MeshSubsets{mesh_items_all_nodes.get()});
        }
        std::vector<unsigned> vec_var_n_components(1, num_components);
        dof_map.reset(
            new NL::LocalToGlobalIndexMap(std::move(components), vec_var_n_components, order));

        auto components_boundary = std::unique_ptr<MeshLib::MeshSubsets>{
            new MeL::MeshSubsets{mesh_items_boundary.get()}};

        dof_map_boundary.reset(dof_map->deriveBoundaryConstrainedMap(
            0,  // variable id
            selected_component,
            std::move(components_boundary),
            boundary_elements));
    }

    template <NL::ComponentOrder order>
    void test(const unsigned num_components, const unsigned selected_component,
              std::function<std::size_t(std::size_t, std::size_t)> const&
                  compute_global_index);

    std::unique_ptr<const MeshLib::Mesh> mesh;
    std::unique_ptr<const MeL::MeshSubset> mesh_items_all_nodes;

    GeoLib::GEOObjects geo_objs;

    std::unique_ptr<NL::LocalToGlobalIndexMap> dof_map;
    std::unique_ptr<NL::LocalToGlobalIndexMap> dof_map_boundary;

    std::unique_ptr<MeL::MeshSubset const> mesh_items_boundary;
    std::vector<MeL::Element*> boundary_elements;
};


struct ComputeGlobalIndexByComponent
{
    std::size_t num_nodes;

    std::size_t operator()(std::size_t const node,
                           std::size_t const component) const
    {
        return node + component * num_nodes;
    }
};

struct ComputeGlobalIndexByLocation
{
    std::size_t num_components;

    std::size_t operator()(std::size_t const node,
                           std::size_t const component) const
    {
        return node * num_components + component;
    }
};


template <NL::ComponentOrder ComponentOrder>
void NumLibLocalToGlobalIndexMapMultiDOFTest::test(
    const unsigned num_components,
    const unsigned selected_component,
    std::function<std::size_t(std::size_t, std::size_t)> const&
        compute_global_index)
{
    initComponents(num_components, selected_component, ComponentOrder);

    ASSERT_EQ(dof_map->getNumberOfComponents(), num_components);
    ASSERT_EQ(dof_map->size(), mesh->getNumberOfElements());

    ASSERT_EQ(dof_map_boundary->getNumberOfComponents(), 1);
    ASSERT_EQ(dof_map_boundary->size(), boundary_elements.size());

    // check mesh elements
    for (unsigned e=0; e<dof_map->size(); ++e)
    {
        auto const element_nodes_size = mesh->getElement(e)->getNumberOfNodes();
        auto const ptr_element_nodes = mesh->getElement(e)->getNodes();

        for (unsigned c=0; c<dof_map->getNumberOfComponents(); ++c)
        {
            auto const& global_idcs = (*dof_map)(e, c).rows;
            ASSERT_EQ(element_nodes_size, global_idcs.size());

            for (unsigned n = 0; n < element_nodes_size; ++n)
            {
                auto const node_id = ptr_element_nodes[n]->getID();
                auto const glob_idx = compute_global_index(node_id, c);
                EXPECT_EQ(glob_idx, global_idcs[n]);
            }
        }
    }

    // check boundary elements
    for (unsigned e=0; e<dof_map_boundary->size(); ++e)
    {
        ASSERT_EQ(1, dof_map_boundary->getNumberOfComponents());

        for (unsigned c=0; c<1; ++c)
        {
            auto const& global_idcs = (*dof_map_boundary)(e, c).rows;

            ASSERT_EQ(2, global_idcs.size()); // boundary of quad is line with two nodes

            for (unsigned n=0; n<2; ++n) // boundary of quad is line with two nodes
            {
                auto const node = e + n;
                auto const glob_idx =
                    compute_global_index(node, selected_component);
                EXPECT_EQ(glob_idx, global_idcs[n]);
            }
        }
    }
}




void assert_equal(NL::LocalToGlobalIndexMap const& dof1, NL::LocalToGlobalIndexMap const& dof2)
{
    ASSERT_EQ(dof1.size(), dof2.size());
    ASSERT_EQ(dof1.getNumberOfComponents(), dof2.getNumberOfComponents());

    for (unsigned e=0; e<dof1.size(); ++e)
    {
        for (unsigned c=0; c<dof1.getNumberOfComponents(); ++c)
        {
            EXPECT_EQ(dof1(e, c).rows, dof2(e, c).rows);
            EXPECT_EQ(dof1(e, c).columns, dof2(e, c).columns);
        }
    }
}

#ifndef USE_PETSC
TEST_F(NumLibLocalToGlobalIndexMapMultiDOFTest, Test1Comp)
#else
TEST_F(NumLibLocalToGlobalIndexMapMultiDOFTest, DISABLED_Test1Comp)
#endif
{
    unsigned const num_components = 1;

    test<NL::ComponentOrder::BY_LOCATION>(
        num_components, 0, ComputeGlobalIndexByComponent{num_components});

    auto dof_map_bc = std::move(dof_map);
    auto dof_map_boundary_bc = std::move(dof_map_boundary);

    test<NL::ComponentOrder::BY_COMPONENT>(
        num_components, 0,
        ComputeGlobalIndexByComponent{(mesh_subdivs + 1) * (mesh_subdivs + 1)});

    assert_equal(*dof_map, *dof_map_bc);
    assert_equal(*dof_map_boundary, *dof_map_boundary_bc);
}

#ifndef USE_PETSC
TEST_F(NumLibLocalToGlobalIndexMapMultiDOFTest, TestMultiCompByComponent)
#else
TEST_F(NumLibLocalToGlobalIndexMapMultiDOFTest, DISABLED_TestMultiCompByComponent)
#endif
{
    unsigned const num_components = 5;
    for (unsigned c = 0; c < num_components; ++c)
        test<NL::ComponentOrder::BY_COMPONENT>(
            num_components, c, ComputeGlobalIndexByComponent{
                                   (mesh_subdivs + 1) * (mesh_subdivs + 1)});
}

#ifndef USE_PETSC
TEST_F(NumLibLocalToGlobalIndexMapMultiDOFTest, TestMultiCompByLocation)
#else
TEST_F(NumLibLocalToGlobalIndexMapMultiDOFTest, DISABLED_TestMultiCompByLocation)
#endif
{
    unsigned const num_components = 5;
    for (unsigned c = 0; c < num_components; ++c)
        test<NL::ComponentOrder::BY_LOCATION>(
            num_components, c, ComputeGlobalIndexByLocation{num_components});
}
