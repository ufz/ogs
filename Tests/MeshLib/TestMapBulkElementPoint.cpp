/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <Eigen/Dense>
#include <array>

#include "MeshLib/Elements/Elements.h"
#include "MeshLib/Elements/MapBulkElementPoint.h"
#include "NumLib/Fem/CoordinatesMapping/NaturalNodeCoordinates.h"
#include "NumLib/Fem/FiniteElement/ElementTraitsLagrange.h"
#include "NumLib/Fem/InitShapeMatrices.h"
#include "Tests/MathLib/PointUtils.h"
#include "Tests/MeshLib/ElementUtils.h"
#include "Tests/NumLib/ReferenceElementUtils.h"
#include "Tests/NumLib/ShapeFunctionUtils.h"

// helpers /////////////////////////////////////////////////////////////////////

namespace
{

// Finds a node in a list of nodes.
template <typename ListOfNodes>
static std::size_t findNode(ListOfNodes const& list_of_nodes,
                            MathLib::Point3d const& node)
{
    auto const it = std::find_if(begin(list_of_nodes),
                                 end(list_of_nodes),
                                 [&node](auto const& coords)
                                 { return MathLib::Point3d{coords} == node; });
    if (it == end(list_of_nodes))
    {
        return std::numeric_limits<std::size_t>::max();
    }

    return std::distance(begin(list_of_nodes), it);
}

// Copies the coordinates of selected nodes in a 3 x #nodes matrix.
//
// natural_coordss is a container of MathLib::Point3d.
// selection is a container of Boolean values, true meaning a node has been
// selected.
template <typename Coords, typename Selection>
static Eigen::MatrixXd getSelectedNodes(Coords const& natural_coordss,
                                        Selection const& selection)
{
    EXPECT_EQ(natural_coordss.size(), selection.size());

    auto const sel_size = std::count(begin(selection), end(selection), true);

    Eigen::MatrixXd sel_coordss(3, sel_size);

    for (std::size_t i_in = 0, i_out = 0; i_in < selection.size(); ++i_in)
    {
        if (!selection[i_in])
        {
            continue;
        }

        Eigen::Map<const Eigen::Vector3d> coords(natural_coordss[i_in].data());

        sel_coordss.col(i_out).noalias() = coords;

        ++i_out;
    }

    return sel_coordss;
}

// Appends the coordinates of a MathLib::Point3d to a 3 x #nodes matrix of node
// coordinates.
static Eigen::MatrixXd appendNode(Eigen::MatrixXd const& mat,
                                  MathLib::Point3d const& pt)
{
    Eigen::Map<const Eigen::Vector3d> col(pt.getCoords());
    Eigen::MatrixXd nodes(mat.rows(), mat.cols() + 1);
    nodes << mat, col;
    return nodes;
}

// Computes the dimension of a flat, i.e., not curved, face from n points on
// that face given as a d_S x n matrix, where d_S is the space dimension.
static unsigned computeFaceDimensionFromPointsOnFace(Eigen::MatrixXd points)
{
    auto const dim = points.rows();
    auto const num_pts = points.cols();
    Eigen::MatrixXd directions_on_face(dim, num_pts - 1);

    for (Eigen::Index i = 0; i < num_pts - 1; ++i)
    {
        directions_on_face.col(i).noalias() = points.col(i + 1) - points.col(0);
    }

    auto const qr_decomposition = directions_on_face.fullPivHouseholderQr();

    // There can only be d_F linearly independent directions on a __flat__ face
    // of dimension d_F. This number d_F is exactly the rank of the directions
    // matrix.
    return static_cast<unsigned>(qr_decomposition.rank());
}
}  // namespace

// the test itself /////////////////////////////////////////////////////////////

template <typename MeshElementType>
class MeshLibMapBulkElementPointTest : public ::testing::Test
{
    ReferenceElementUtils::ReferenceElement<MeshElementType> reference_element;

protected:
    MeshElementType const& bulk_element = reference_element.element;
};

// The mesh element types for which we test the coordinate mapping.
using MeshElementTypes = ::testing::Types<MeshLib::Line,
                                          MeshLib::Quad,
                                          MeshLib::Hex,
                                          MeshLib::Tri,
                                          MeshLib::Tet,
                                          MeshLib::Prism,
                                          MeshLib::Pyramid>;

TYPED_TEST_SUITE(MeshLibMapBulkElementPointTest, MeshElementTypes);

TYPED_TEST(MeshLibMapBulkElementPointTest, Prerequisites)
{
    auto const num_faces = ElementUtils::getNumberOfFaces(this->bulk_element);

    // TODO this is a test of ElementUtils;
    ASSERT_LT(0, num_faces)
        << "Each tested element type should have at least one face. Note, in "
           "this test suite 'faces' denote the d-1 dimensional faces of the d "
           "dimensional bulk element.";
}

template <typename MeshElementType>
static void mapAllFaceNodesToBulkElementAndCheckIfTheyAreFound(
    std::size_t const face_id,
    MeshElementType const& bulk_element,
    std::array<bool, MeshElementType::n_all_nodes>&
        has_bulk_element_node_been_visited,
    bool& some_node_of_this_face_has_not_been_found_in_the_bulk_element)
{
    has_bulk_element_node_been_visited =
        std::array<bool, MeshElementType::n_all_nodes>{};
    some_node_of_this_face_has_not_been_found_in_the_bulk_element = false;

    auto const bulk_element_node_natural_coords =
        ReferenceElementUtils::getNodeCoordsOfReferenceElement<
            MeshElementType>();

    auto const face_ptr = ElementUtils::getFace(bulk_element, face_id);
    auto const face_cell_type = face_ptr->getCellType();
    auto const face_nodes_natural_coordss =
        ReferenceElementUtils::getNodeCoordsOfReferenceElement(face_cell_type);

    for (auto const& face_node_natural_coords : face_nodes_natural_coordss)
    {
        auto const wp = PointUtils::getWeightedPointOnFace<MeshElementType>(
            face_node_natural_coords);

        ASSERT_EQ(face_ptr->getDimension(), wp.getDimension())
            << "Unexpected dimension. This is probably an internal logic "
               "error of this unit test. The right dimension is crucial, "
               "otherwise the wrong mapping function would be used.";

        // This is the actual testee.
        auto const face_node_mapped_to_bulk_element =
            MeshLib::getBulkElementPoint(bulk_element, face_id, wp);

        // We try to find the mapped natural coordinates in the list of bulk
        // element corners. This assumes that natural coordinates coincide with
        // real coordinates, which holds true for many element types, but not
        // for pyramids.
        auto const found_bulk_node_index = findNode(
            bulk_element_node_natural_coords, face_node_mapped_to_bulk_element);

        if (found_bulk_node_index == std::numeric_limits<std::size_t>::max())
        {
            ADD_FAILURE() << "A mapped node has not been found in the list "
                             "of bulk element nodes."
                             "\n  face id:          "
                          << face_id
                          << "\n  face node:        " << std::setprecision(16)
                          << MathLib::Point3d(face_node_natural_coords)
                          << " (dim: " << face_ptr->getDimension() << ')'
                          << "\n  mapped bulk node: "
                          << face_node_mapped_to_bulk_element
                          << " (not a node of the bulk element)";
            some_node_of_this_face_has_not_been_found_in_the_bulk_element =
                true;
            continue;
        }

        EXPECT_FALSE(has_bulk_element_node_been_visited[found_bulk_node_index])
            << "Same node has been found twice for a single face. That "
               "must not happen, because different face node must be "
               "mapped to different bulk element nodes.";

        has_bulk_element_node_been_visited[found_bulk_node_index] = true;
    }
}

// This test checks consistency between the face-to-bulk natural coordinates
// mapping and the natural coordinates used in FEM shape functions.
//
// We iterate over all nodes of all faces of a bulk element and check that these
// nodes are mapped to nodes of the bulk element. I.e., this test checks if the
// face nodes are mapped to some nodes of the bulk element; it does not check if
// the face nodes are mapped to the right bulk nodes.
TYPED_TEST(MeshLibMapBulkElementPointTest,
           CheckFaceNodesAreMappedToBulkElementNodes)
{
    using MeshElementType = TypeParam;

    if constexpr (std::is_same_v<MeshElementType, MeshLib::Pyramid>)
    {
        // This test assumes that for reference elements the natural coordinates
        // coincide with the real coordinates. This is true for most of OGS's
        // element types, but not for pyramids, whose natural coordinates are
        // from [-1,1]^3 but whose volume covers only a subset of that cube.
        return;  // skip this test quietly
    }

    auto const bulk_element_node_natural_coords =
        ReferenceElementUtils::getNodeCoordsOfReferenceElement<
            MeshElementType>();
    auto const bulk_element_num_nodes = MeshElementType::n_all_nodes;

    auto const num_faces = ElementUtils::getNumberOfFaces(this->bulk_element);

    // Saves how many times face nodes have been mapped to a specific bulk
    // element node.
    std::array<unsigned, bulk_element_num_nodes> num_node_visits_all_faces{};

    bool some_node_of_any_face_has_not_been_found_in_the_bulk_element =
        false;  // for diagnostic information only

    for (std::size_t face_id = 0; face_id < num_faces; ++face_id)
    {
        // Saves if a specific bulk element node has been mapped to from some
        // node of the current face.
        std::array<bool, bulk_element_num_nodes> has_node_been_visited{};
        bool some_node_of_this_face_has_not_been_found_in_the_bulk_element =
            false;

        mapAllFaceNodesToBulkElementAndCheckIfTheyAreFound(
            face_id,
            this->bulk_element,
            has_node_been_visited,
            some_node_of_this_face_has_not_been_found_in_the_bulk_element);

        // Collect data of all faces
        if (some_node_of_this_face_has_not_been_found_in_the_bulk_element)
        {
            some_node_of_any_face_has_not_been_found_in_the_bulk_element = true;
        }

        for (unsigned n = 0; n < bulk_element_num_nodes; ++n)
        {
            if (has_node_been_visited[n])
            {
                ++num_node_visits_all_faces[n];
            }
        }
    }

    // After iterating over all nodes of all faces we should have
    // found/visited/mapped each node of the bulk element at least once.
    // Note: This is not necessarily true for higher order mesh elements.
    for (std::size_t i = 0; i < bulk_element_num_nodes; ++i)
    {
        auto const num_visits = num_node_visits_all_faces[i];
        EXPECT_NE(0, num_visits)
            << "Bulk node #" << i << " has never been found for any face.";
    }

    // Print some diagnostic information
    if (some_node_of_any_face_has_not_been_found_in_the_bulk_element)
    {
        ADD_FAILURE()
            << "Some face node has not been mapped to a bulk element node.";

        std::cout << "Bulk element's nodes' natural coordinates are:\n";
        for (std::size_t i = 0; i < bulk_element_num_nodes; ++i)
        {
            std::cout << "  #" << i << ":\t"
                      << MathLib::Point3d(bulk_element_node_natural_coords[i])
                      << '\n';
        }
    }
}

template <typename MeshElementType>
static void checkInterpolatedFaceNodes(
    Eigen::MatrixXd const& visited_bulk_element_nodes,
    std::size_t const face_id, MeshElementType const& bulk_element)
{
    auto constexpr face_dimension = MeshElementType::dimension - 1;

    auto const face_dimension_from_points =
        computeFaceDimensionFromPointsOnFace(visited_bulk_element_nodes);

    // Check internal consistency - dimensionality of the current face.
    // Note: this check does not work for elements with curved faces.
    ASSERT_EQ(face_dimension, face_dimension_from_points)
        << "The face dimension computed from the face's nodes does not match "
           "the dimension obtained from the face's mesh element type. This "
           "might be due to floating point accuracy or due to internal logic "
           "errors of this unit test.";

    auto const face_ptr = ElementUtils::getFace(bulk_element, face_id);

    // Add an additional point on the face to the list of bulk element nodes and
    // recompute the dimension of the face.
    for (auto const& natural_coords :
         ReferenceElementUtils::getCoordsInReferenceElementForTest(*face_ptr))
    {
        auto const wp =
            PointUtils::getWeightedPointOnFace<MeshElementType>(natural_coords);

        // This is the actual testee.
        auto const face_node_mapped_to_bulk_element =
            MeshLib::getBulkElementPoint(bulk_element, face_id, wp);

        // This check assumes that natural coordinates and real coordinates of
        // the bulk element are related. I.e., either they coincide or the
        // natural coordinates are from a subset of the real coordinates domain.
        // Both is true for most of OGS's element types, none is true for
        // pyramids.
        ASSERT_TRUE(bulk_element.isPntInElement(
            face_node_mapped_to_bulk_element, 1e-14))
            << "Mapped point is not in bulk element. This could also be an "
               "issue with floating point precision."
               "\nwp on face "
            << wp  //
            << "\nmapped pt  " << face_node_mapped_to_bulk_element;

        auto const nodes_natural_coords = appendNode(
            visited_bulk_element_nodes, face_node_mapped_to_bulk_element);

        auto const face_dimension_from_points2 =
            computeFaceDimensionFromPointsOnFace(nodes_natural_coords);

        // Note: this check does not work for elements with curved faces.
        EXPECT_EQ(face_dimension, face_dimension_from_points2)
            << "Either the added point does not lie on the face. Or there are "
               "issues with numerical precision or with the internal logic of "
               "this test case."
               "\nwp on face "
            << wp  //
            << "\nmapped pt  " << face_node_mapped_to_bulk_element;
    }
}

// This test is very similar to CheckFaceNodesAreMappedToBulkElementNodes above.
//
// It checks a subset of the test above and extends it with a check that points
// that are not corners of the face are mapped to a bulk point lying on the
// surface of the bulk element.
TYPED_TEST(MeshLibMapBulkElementPointTest, CheckInterpolatedFaceNodes)
{
    using MeshElementType = TypeParam;

    if constexpr (std::is_same_v<MeshElementType, MeshLib::Pyramid>)
    {
        // This test requires that the coordinates of all bulk element nodes are
        // in the range (i.e., result set) of MeshLib::getBulkElementPoint().
        // That implicitly means that natural coordinates and real coordinates
        // must coincide for reference elements. That, in turn, is not true for
        // pyramids.
        return;  // skip test quietly
    }

    auto const bulk_element_node_natural_coords =
        ReferenceElementUtils::getNodeCoordsOfReferenceElement<
            MeshElementType>();
    auto const bulk_element_num_nodes = MeshElementType::n_all_nodes;

    auto const num_faces = ElementUtils::getNumberOfFaces(this->bulk_element);

    for (std::size_t face_id = 0; face_id < num_faces; ++face_id)
    {
        // Saves if a specific bulk element node has been mapped to from some
        // node of the current face.
        std::array<bool, bulk_element_num_nodes> has_node_been_visited{};
        bool some_node_of_this_face_has_not_been_found_in_the_bulk_element =
            false;

        mapAllFaceNodesToBulkElementAndCheckIfTheyAreFound(
            face_id,
            this->bulk_element,
            has_node_been_visited,
            some_node_of_this_face_has_not_been_found_in_the_bulk_element);

        if (some_node_of_this_face_has_not_been_found_in_the_bulk_element)
        {
            // The remainder of the loop body only makes sense if the part above
            // was successful.
            continue;
        }

        auto const visited_bulk_element_nodes = getSelectedNodes(
            bulk_element_node_natural_coords, has_node_been_visited);

        checkInterpolatedFaceNodes(visited_bulk_element_nodes, face_id,
                                   this->bulk_element);
    }
}

// Maps and checks all corners of a face to the bulk element.
template <typename BulkElementType>
static void mapAllFaceNodesToBulkElementAndCheckWithShapeFunctionInterpolation(
    std::size_t const face_id, BulkElementType const& bulk_element)
{
    auto const bulk_element_node_natural_coords =
        ReferenceElementUtils::getNodeCoordsOfReferenceElement<
            BulkElementType>();
    auto const bulk_nodes_wps =
        PointUtils::toWeightedPointsOfDim<BulkElementType::dimension>(
            bulk_element_node_natural_coords);
    auto const shape_matrices_at_bulk_nodes =
        ShapeFunctionUtils::computeShapeMatricesBulk(bulk_element,
                                                     bulk_nodes_wps);

    auto const face_ptr = ElementUtils::getFace(bulk_element, face_id);
    auto const face_cell_type = face_ptr->getCellType();
    auto const face_nodes_natural_coordss =
        ReferenceElementUtils::getNodeCoordsOfReferenceElement(face_cell_type);
    auto const num_face_nodes = face_ptr->getNumberOfNodes();

    auto const face_nodes_wps =
        PointUtils::toWeightedPointsOnFace<BulkElementType>(
            face_nodes_natural_coordss);
    auto const shape_matrices_at_face_nodes =
        ShapeFunctionUtils::computeShapeMatricesOnFace(*face_ptr,
                                                       face_nodes_wps);

    for (std::size_t i_face_node = 0; i_face_node < num_face_nodes;
         ++i_face_node)
    {
        auto const i_bulk_node = ElementUtils::getLocalBulkNodeIndexOfFaceNode(
            bulk_element, face_id, i_face_node);
        auto const& wp_in_face = face_nodes_wps[i_face_node];

        // This is the actual testee.
        auto const natural_coords_face_node_in_bulk_element =
            MeshLib::getBulkElementPoint(bulk_element, face_id, wp_in_face);

        auto const& N_face = shape_matrices_at_face_nodes[i_face_node].N;

        // This is the reference value. It is obtained without relying on the
        // mapping procedure under test (of course).
        auto const node_coords_from_face_interpolation =
            ShapeFunctionUtils::interpolateNodeCoordinatesFace(*face_ptr,
                                                               N_face);

        auto const& N_bulk = shape_matrices_at_bulk_nodes[i_bulk_node].N;
        auto const node_coords_from_bulk_interpolation =
            ShapeFunctionUtils::interpolateNodeCoordinates(bulk_element,
                                                           N_bulk);

        // This is __the__ fundamental property of the face to bulk mapping.
        EXPECT_EQ(node_coords_from_face_interpolation,
                  node_coords_from_bulk_interpolation)
            << "Interpolation on the face must yield the same result as "
               "interpolation in the bulk element."
            << "\nface         #" << face_id      //
            << "\nnode of face #" << i_face_node  //
            << "\nnode of bulk #" << i_bulk_node;
    }
}

// Maps and checks a variety of points on a face to the bulk element.
template <typename BulkElementType>
static void mapPointsOnFaceToBulkElementAndCheckWithShapeFunctionInterpolation(
    std::size_t const face_id, BulkElementType const& bulk_element)
{
    auto const face_ptr = ElementUtils::getFace(bulk_element, face_id);
    auto const natural_coords_face =
        ReferenceElementUtils::getCoordsInReferenceElementForTest(*face_ptr);
    auto const wps_face = PointUtils::toWeightedPointsOnFace<BulkElementType>(
        natural_coords_face);
    auto const shape_matrices_face =
        ShapeFunctionUtils::computeShapeMatricesOnFace(*face_ptr, wps_face);

    for (std::size_t i = 0; i < shape_matrices_face.size(); ++i)
    {
        auto const& wp_in_face = wps_face[i];
        auto const& N_face = shape_matrices_face[i].N;

        // This is the reference value. It is obtained without relying on the
        // mapping procedure under test (of course).
        auto const node_coords_from_face_interpolation =
            ShapeFunctionUtils::interpolateNodeCoordinatesFace(*face_ptr,
                                                               N_face);

        // This is the actual testee.
        auto const face_point_in_bulk_element =
            MeshLib::getBulkElementPoint(bulk_element, face_id, wp_in_face);

        auto const bulk_natural_coords =
            PointUtils::getCoords(face_point_in_bulk_element);

        auto const wps_in_bulk =
            PointUtils::toWeightedPointsOfDim<BulkElementType::dimension>(
                std::array<std::array<double, 3>, 1>{bulk_natural_coords});

        auto const N_bulk = ShapeFunctionUtils::computeShapeMatricesBulk(
                                bulk_element, wps_in_bulk)
                                .front()
                                .N;

        auto const node_coords_from_bulk_interpolation =
            ShapeFunctionUtils::interpolateNodeCoordinates(bulk_element,
                                                           N_bulk);

        // This is __the__ fundamental property of the face to bulk mapping.
        ASSERT_PRED_FORMAT3(PointUtils::IsNear{},
                            node_coords_from_face_interpolation,
                            node_coords_from_bulk_interpolation, 1e-15)
            << "Interpolation of the face nodes on the face must yield the same"
               "result as interpolation of the bulk nodes in the bulk element."
            << "\nface         #" << face_id     //
            << "\nwp in face    " << wp_in_face  //
            << "\nwp in bulk    " << wps_in_bulk.front();
    }
}

// This test is much stricter than the tests before. It does not only check if
// the face to bulk mapping under test maps face nodes to some bulk element
// nodes, but also if they are mapped to the right bulk element nodes. However,
// this test involves shape function interpolation, i.e., needs a larger
// machinery. Therefore, the other tests are still justified.
TYPED_TEST(MeshLibMapBulkElementPointTest,
           CheckFaceNodesAreMappedToTheRightBulkElementNodes)
{
    auto const num_faces = ElementUtils::getNumberOfFaces(this->bulk_element);

    for (std::size_t face_id = 0; face_id < num_faces; ++face_id)
    {
        mapAllFaceNodesToBulkElementAndCheckWithShapeFunctionInterpolation(
            face_id, this->bulk_element);
    }
}

// Same for points that are anywhere on the face, not only (corner) nodes of the
// face.
TYPED_TEST(MeshLibMapBulkElementPointTest,
           CheckPointsOnFaceCoincideWithMappedPoints)
{
    auto const num_faces = ElementUtils::getNumberOfFaces(this->bulk_element);

    for (std::size_t face_id = 0; face_id < num_faces; ++face_id)
    {
        mapPointsOnFaceToBulkElementAndCheckWithShapeFunctionInterpolation(
            face_id, this->bulk_element);
    }
}

TEST(MeshLib, MapBulkElementPoint3D)
{
    MathLib::WeightedPoint wp{std::array{0.1, 0.2, 0.3}, 0.4};
    auto const mapped_point =
        MeshLib::getBulkElementPoint(MeshLib::CellType::HEX8, 0, wp);

    // The mapping for "3D faces" is special as it will always return the
    // origin.
    ASSERT_EQ(MathLib::ORIGIN, mapped_point);
}
