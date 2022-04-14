/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "Tests/MathLib/PointUtils.h"
#include "Tests/NumLib/ReferenceElementUtils.h"
#include "Tests/NumLib/ShapeFunctionUtils.h"

template <typename MeshElementType>
class NumLibReferenceElementTest : public ::testing::Test
{
    ReferenceElementUtils::ReferenceElement<MeshElementType> reference_element;

protected:
    MeshElementType const& element = reference_element.element;
};

using MeshElementTypes =
    ConvertListType_t<MeshLib::AllElementTypes, ::testing::Types>;

TYPED_TEST_SUITE(NumLibReferenceElementTest, MeshElementTypes);

template <typename MeshElementType, typename NaturalCoordsContainer>
static void interpolateNodeCoordsAndCheckTheResult(
    MeshElementType const& element,
    NaturalCoordsContainer const& natural_coords)
{
    auto const wps =
        PointUtils::toWeightedPointsOfDim<MeshElementType::dimension>(
            natural_coords);
    auto const shape_matrices =
        ShapeFunctionUtils::computeShapeMatricesBulk(element, wps);

    for (std::size_t i = 0; i < shape_matrices.size(); ++i)
    {
        auto const& N = shape_matrices[i].N;

        auto const node_coords_from_bulk_interpolation =
            ShapeFunctionUtils::interpolateNodeCoordinates(element, N);

        // Only holds if the natural coordinates coincide with the real
        // coordinates on the reference element. For pyramids this is true only
        // for certain coordinates.
        EXPECT_TRUE(
            element.isPntInElement(node_coords_from_bulk_interpolation, 1e-14))
            << "The interpolated node is not inside the element. This "
               "could "
               "also be an issue with numerical precision."
               "\n i = "
            << i;

        // Only holds if the natural coordinates coincide with the real
        // coordinates on the reference element. For pyramids this is true only
        // for certain coordinates.
        ASSERT_PRED_FORMAT3(PointUtils::IsNear{}, natural_coords[i],
                            node_coords_from_bulk_interpolation, 1e-15)
            << "The interpolation takes place on the reference element, so the "
               "coordinates interpolated from the bulk element nodes must be "
               "equal to the natural coordinates."
               "\ni = "
            << i;
    }
}

// TODO this is a test of ReferenceElementUtils
TYPED_TEST(NumLibReferenceElementTest, Prerequisite)
{
    using MeshElementType = TypeParam;

    auto const element_node_natural_coords =
        ReferenceElementUtils::getNodeCoordsOfReferenceElement<
            MeshElementType>();

    ASSERT_EQ(MeshElementType::n_all_nodes, element_node_natural_coords.size())
        << "Number of nodes and natural coordinates do not match for this "
           "element. This could be an internal logic error of this unit test.";
}

// For most element types the element's natural coordinates coincide with the
// real coordinates in the reference element. I.e., when interpolating the node
// coordinates of the reference element to a point inside the element given by
// its natural coordinates r, s, t the following holds:
//
// (r, s, t) == N(r, s, t) * [coords of reference element's nodes],
//
// where N is a combination of shape matrices used for interpolation.
//
// One notable exception are pyramids.
//
// This unit test checks the property described above for r, s, t being the
// coordinates of the nodes of the reference element.
//
TYPED_TEST(NumLibReferenceElementTest, ElementNodes)
{
    using MeshElementType = TypeParam;

    auto const all_element_node_natural_coords =
        ReferenceElementUtils::getNodeCoordsOfReferenceElement<
            MeshElementType>();

    // For pyramid 13 some higher order nodes would not pass this unit test.
    auto const element_node_natural_coords =
        !std::is_same_v<MeshElementType, MeshLib::Pyramid13>
            ? all_element_node_natural_coords
            : BaseLib::DynamicSpan{all_element_node_natural_coords.data,
                                   MeshLib::Pyramid13::n_base_nodes};

    interpolateNodeCoordsAndCheckTheResult(this->element,
                                           element_node_natural_coords);
}

// For most element types the element's natural coordinates coincide with the
// real coordinates in the reference element. I.e., when interpolating the node
// coordinates of the reference element to a point inside the element given by
// its natural coordinates r, s, t the following holds:
//
// (r, s, t) == N(r, s, t) * [coords of reference element's nodes],
//
// where N is a combination of shape matrices used for interpolation.
//
// One notable exception are pyramids.
//
// This unit test checks the property described above for r, s, t being the
// coordinates some points inside the reference element.
//
TYPED_TEST(NumLibReferenceElementTest, PointsInsideElement)
{
    using MeshElementType = TypeParam;

    if constexpr (std::is_same_v<MeshElementType, MeshLib::Pyramid> ||
                  std::is_same_v<MeshElementType, MeshLib::Pyramid13>)
    {
        GTEST_SKIP() << "Test does not apply to pyramids.";
    }

    auto const natural_coords =
        ReferenceElementUtils::getCoordsInReferenceElementForTest(
            this->element);

    interpolateNodeCoordsAndCheckTheResult(this->element, natural_coords);
}
