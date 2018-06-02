/*
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <gtest/gtest.h>

#include <memory>

#include "GeoLib/GEOObjects.h"
#include "GeoLib/Polyline.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/Location.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/MeshGenerators/MeshGenerator.h"
#include "MeshLib/MeshSearch/NodeSearch.h"
#include "MeshLib/MeshSubset.h"
#include "NumLib/DOF/LocalToGlobalIndexMap.h"

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
        mesh_items_all_nodes =
            std::make_unique<MeL::MeshSubset>(*mesh, mesh->getNodes());

        auto ply_pnts = std::make_unique<std::vector<GeoLib::Point*>>();
        ply_pnts->push_back(new GeoLib::Point(0.0, 0.0, 0.0));
        ply_pnts->push_back(new GeoLib::Point(1.0, 0.0, 0.0));

        std::string geometry_0("GeometryWithPntsAndPolyline");
        geo_objs.addPointVec(std::move(ply_pnts), geometry_0, nullptr);

        auto ply = new GeoLib::Polyline(*geo_objs.getPointVec(geometry_0));

        ply->addPoint(0);
        ply->addPoint(1);

        auto plys = std::make_unique<std::vector<GeoLib::Polyline*>>();
        plys->push_back(ply);

        geo_objs.addPolylineVec(std::move(plys), geometry_0, nullptr);

        auto search_length = std::make_unique<MGTL::SearchLength>();
        MGTL::MeshNodeSearcher const& searcher_nodes =
            MGTL::MeshNodeSearcher::getMeshNodeSearcher(
                *mesh, std::move(search_length));
        MGTL::BoundaryElementsSearcher searcher_elements(*mesh, searcher_nodes);

        auto elems = searcher_elements.getBoundaryElements(*ply);

        // deep copy because the searcher destroys the elements.
        std::transform(elems.cbegin(), elems.cend(),
                       std::back_inserter(boundary_elements),
                       [](MeL::Element const* e) { return e->clone(); });

        std::vector<MeL::Node*> nodes = MeL::getUniqueNodes(boundary_elements);

        nodes_subset =
            nodesNodesIntersection(mesh_items_all_nodes->getNodes(), nodes);
        mesh_items_boundary = std::make_unique<MeshLib::MeshSubset>(
            mesh_items_all_nodes->getMesh(), nodes_subset);
    }

    ~NumLibLocalToGlobalIndexMapMultiDOFTest() override
    {
        for (auto e : boundary_elements)
            delete e;
    }

    void initComponents(const int num_components,
                        const int selected_component,
                        const NL::ComponentOrder order)
    {
        assert(selected_component < static_cast<int>(num_components));

        std::vector<MeshLib::MeshSubset> components;
        for (int i = 0; i < num_components; ++i)
        {
            components.emplace_back(*mesh_items_all_nodes);
        }
        std::vector<int> vec_var_n_components(1, num_components);
        dof_map = std::make_unique<NL::LocalToGlobalIndexMap>(
            std::move(components), vec_var_n_components, order);

        MeL::MeshSubset components_boundary{*mesh_items_boundary};

        dof_map_boundary.reset(dof_map->deriveBoundaryConstrainedMap(
            0,  // variable id
            {selected_component},
            std::move(components_boundary),
            boundary_elements));
    }

    // Multi-component version.
    void initComponents(int const num_components,
                        std::vector<int> const& selected_components,
                        NL::ComponentOrder const order)
    {
        assert(static_cast<int>(selected_components.size()) <= num_components);

        std::vector<MeshLib::MeshSubset> components;
        for (int i = 0; i < num_components; ++i)
        {
            components.emplace_back(*mesh_items_all_nodes);
        }
        std::vector<int> vec_var_n_components(1, num_components);
        dof_map = std::make_unique<NL::LocalToGlobalIndexMap>(
            std::move(components), vec_var_n_components, order);

        MeL::MeshSubset components_boundary{*mesh_items_boundary};

        dof_map_boundary.reset(dof_map->deriveBoundaryConstrainedMap(
            0,  // variable id
            selected_components,
            std::move(components_boundary),
            boundary_elements));
    }

    template <NL::ComponentOrder order>
    void test(const int num_components, const int selected_component,
              std::function<std::size_t(std::size_t, std::size_t)> const&
                  compute_global_index);

    // Multicomponent version
    template <NL::ComponentOrder order>
    void test(const int num_components,
              std::vector<int> const& selected_component,
              std::function<std::size_t(std::size_t, std::size_t)> const&
                  compute_global_index);

    std::unique_ptr<const MeshLib::Mesh> mesh;
    std::unique_ptr<MeL::MeshSubset const> mesh_items_all_nodes;

    GeoLib::GEOObjects geo_objs;

    std::unique_ptr<NL::LocalToGlobalIndexMap> dof_map;
    std::unique_ptr<NL::LocalToGlobalIndexMap> dof_map_boundary;

    std::unique_ptr<MeL::MeshSubset const> mesh_items_boundary;
    std::vector<MeL::Element*> boundary_elements;

    /// Intersection of boundary nodes and bulk mesh subset.
    std::vector<MeshLib::Node*> nodes_subset;
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
    const int num_components,
    const int selected_component,
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

        for (int c = 0; c < dof_map->getNumberOfComponents(); ++c)
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

        for (int c = 0; c < 1; ++c)
        {
            auto const& global_idcs = (*dof_map_boundary)(e, c).rows;

            ASSERT_EQ(2, global_idcs.size()); // boundary of quad is line with two nodes

            // e is the number of a boundary element (which is a line, so two
            // nodes) from 0 to something. e+n must be the node ids along the
            // x-axis of the mesh;
            for (unsigned n = 0; n < 2; ++n)
            {
                auto const node = e + n;
                auto const glob_idx =
                    compute_global_index(node, selected_component);
                EXPECT_EQ(glob_idx, global_idcs[n]);
            }
        }
    }
}

// Multicomponent version
template <NL::ComponentOrder ComponentOrder>
void NumLibLocalToGlobalIndexMapMultiDOFTest::test(
    int const num_components,
    std::vector<int> const& selected_components,
    std::function<std::size_t(std::size_t, std::size_t)> const&
        compute_global_index)
{
    initComponents(num_components, selected_components, ComponentOrder);

    ASSERT_EQ(dof_map->getNumberOfComponents(), num_components);
    ASSERT_EQ(dof_map->size(), mesh->getNumberOfElements());

    ASSERT_EQ(dof_map_boundary->getNumberOfComponents(),
              selected_components.size());
    ASSERT_EQ(dof_map_boundary->size(), boundary_elements.size());

    // check mesh elements
    for (unsigned e = 0; e < dof_map->size(); ++e)
    {
        auto const element_nodes_size = mesh->getElement(e)->getNumberOfNodes();
        auto const ptr_element_nodes = mesh->getElement(e)->getNodes();

        for (int c = 0; c < dof_map->getNumberOfComponents(); ++c)
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
    for (unsigned e = 0; e < dof_map_boundary->size(); ++e)
    {
        ASSERT_EQ(selected_components.size(),
                  dof_map_boundary->getNumberOfComponents());

        for (int c = 0; c < static_cast<int>(selected_components.size()); ++c)
        {
            auto const& global_idcs = (*dof_map_boundary)(e, c).rows;

            ASSERT_EQ(
                2,
                global_idcs.size());  // boundary of quad is line with two nodes

            // e is the number of a boundary element (which is a line, so two
            // nodes) from 0 to something. e+n must be the node ids along the
            // x-axis of the mesh;
            for (unsigned n = 0; n < 2; ++n)
            {
                auto const node = e + n;
                auto const glob_idx =
                    compute_global_index(node, selected_components[c]);
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
        for (int c = 0; c < dof1.getNumberOfComponents(); ++c)
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
    int const num_components = 1;

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
    int const num_components = 5;
    for (int c = 0; c < num_components; ++c)
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
    int const num_components = 5;
    for (int c = 0; c < num_components; ++c)
        test<NL::ComponentOrder::BY_LOCATION>(
            num_components, c, ComputeGlobalIndexByLocation{num_components});
}

#ifndef USE_PETSC
TEST_F(NumLibLocalToGlobalIndexMapMultiDOFTest, TestMultiComp2ByComponent)
#else
TEST_F(NumLibLocalToGlobalIndexMapMultiDOFTest,
       DISABLED_TestMultiComp2ByComponent)
#endif
{
    int const num_components = 5;
    for (int c = 0; c < num_components; ++c)
    {
        std::vector<int> const components = {c, (c + 1) % num_components};
        test<NL::ComponentOrder::BY_COMPONENT>(
            num_components, components,
            ComputeGlobalIndexByComponent{(mesh_subdivs + 1) *
                                          (mesh_subdivs + 1)});
    }
}
#ifndef USE_PETSC
TEST_F(NumLibLocalToGlobalIndexMapMultiDOFTest, TestMultiComp2ByLocation)
#else
TEST_F(NumLibLocalToGlobalIndexMapMultiDOFTest,
       DISABLED_TestMultiComp2ByLocation)
#endif
{
    int const num_components = 5;
    for (int c = 0; c < num_components; ++c)
    {
        std::vector<int> const components = {c, (c + 1) % num_components};
        test<NL::ComponentOrder::BY_LOCATION>(
            num_components, components,
            ComputeGlobalIndexByLocation{num_components});
    }
}
