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
#include "MeshLib/Elements/Utils.h"
#include "Tests/MeshLib/ElementUtils.h"
#include "Tests/NumLib/ReferenceElementUtils.h"

template <typename MeshElementType>
class MeshLibSurfaceNormalsTest : public ::testing::Test
{
    ReferenceElementUtils::ReferenceElement<MeshElementType> reference_element;

protected:
    MeshElementType const& bulk_element = reference_element.element;
};

using MeshElementTypes =
    ConvertListType_t<MeshLib::AllElementTypes, ::testing::Types>;

TYPED_TEST_SUITE(MeshLibSurfaceNormalsTest, MeshElementTypes);

TYPED_TEST(MeshLibSurfaceNormalsTest, RightDirection)
{
    auto const& e = this->bulk_element;
    auto const n_faces = ElementUtils::getNumberOfFaces(e);

    for (std::size_t face_id = 0; face_id < n_faces; ++face_id)
    {
        auto const face_ptr = ElementUtils::getFace(e, face_id);

        auto const outward_surface_normal_vec =
            calculateNormalizedSurfaceNormal(*face_ptr, e);

        auto const face_center = getCenterOfGravity(*face_ptr);

        Eigen::Vector3d const face_center_vec(face_center[0], face_center[1],
                                              face_center[2]);

        auto const eps = 1e-3;

        {
            Eigen::Vector3d const internal_point_vec =
                face_center_vec - eps * outward_surface_normal_vec;

            MathLib::Point3d const internal_point({internal_point_vec[0],
                                                   internal_point_vec[1],
                                                   internal_point_vec[2]});

            if (!e.isPntInElement(internal_point))
            {
                ADD_FAILURE() << "The computed internal point is outside the "
                                 "bulk element."
                                 "\nThe face id is "
                              << face_id << "\nThe face is " << *face_ptr
                              << "\nThe internal point is " << internal_point
                              << "\nThe outward surface normal is "
                              << outward_surface_normal_vec.transpose();

                if (e.getDimension() == 2)
                {
                    std::cout
                        << "The bulk element normal is "
                        << MeshLib::FaceRule::getSurfaceNormal(e).transpose()
                        << '\n';
                }
            }
        }

        {
            Eigen::Vector3d const external_point_vec =
                face_center_vec + eps * outward_surface_normal_vec;

            MathLib::Point3d const external_point({external_point_vec[0],
                                                   external_point_vec[1],
                                                   external_point_vec[2]});

            if (e.isPntInElement(external_point))
            {
                ADD_FAILURE()
                    << "The computed external point is inside the bulk element."
                       "\nThe face id is "
                    << face_id << "\nThe face is " << *face_ptr
                    << "\nThe external point is " << external_point
                    << "\nThe outward surface normal is "
                    << outward_surface_normal_vec.transpose();

                if (e.getDimension() == 2)
                {
                    std::cout
                        << "The bulk element normal is "
                        << MeshLib::FaceRule::getSurfaceNormal(e).transpose()
                        << '\n';
                }
            }
        }
    }
}

TYPED_TEST(MeshLibSurfaceNormalsTest, UnitLength)
{
    auto const& e = this->bulk_element;
    auto const n_faces = ElementUtils::getNumberOfFaces(e);

    for (std::size_t face_id = 0; face_id < n_faces; ++face_id)
    {
        auto const face_ptr = ElementUtils::getFace(e, face_id);

        auto const outward_surface_normal_vec =
            calculateNormalizedSurfaceNormal(*face_ptr, e);

        EXPECT_DOUBLE_EQ(1.0, outward_surface_normal_vec.norm());
    }
}

TYPED_TEST(MeshLibSurfaceNormalsTest, PerpendicularToSurface)
{
    using MeshElementType = TypeParam;
    constexpr auto dim = MeshElementType::dimension;

    auto const& e = this->bulk_element;
    [[maybe_unused]] auto const n_faces = ElementUtils::getNumberOfFaces(e);

    if constexpr (dim == 0)
    {
        return;
    }
    else if constexpr (dim == 1)
    {
        auto const& a = e.getNode(0)->asEigenVector3d();
        auto const& b = e.getNode(1)->asEigenVector3d();
        Eigen::Vector3d const v = b - a;

        for (std::size_t face_id = 0; face_id < n_faces; ++face_id)
        {
            auto const face_ptr = ElementUtils::getFace(e, face_id);

            auto const outward_surface_normal_vec =
                calculateNormalizedSurfaceNormal(*face_ptr, e);

            double const v_dot_n = v.dot(outward_surface_normal_vec);
            EXPECT_DOUBLE_EQ(v.squaredNorm(), v_dot_n * v_dot_n)
                << "The normal is not parallel to the line element."
                   "\nThe face id is "
                << face_id << "\nThe line's direction is " << v.transpose()
                << "\nThe surface normal is "
                << outward_surface_normal_vec.transpose();
        }
    }
    else if constexpr (dim == 2)
    {
        auto const bulk_normal = MeshLib::FaceRule::getSurfaceNormal(e);

        for (std::size_t face_id = 0; face_id < n_faces; ++face_id)
        {
            auto const face_ptr = ElementUtils::getFace(e, face_id);

            auto const outward_surface_normal_vec =
                calculateNormalizedSurfaceNormal(*face_ptr, e);

            auto const& a = face_ptr->getNode(0)->asEigenVector3d();
            auto const& b = face_ptr->getNode(1)->asEigenVector3d();
            Eigen::Vector3d const v = b - a;

            EXPECT_DOUBLE_EQ(0, v.dot(outward_surface_normal_vec))
                << "The surface normal is not perpendicular to the edge."
                   "\nThe face id is "
                << face_id << "\nThe edge's direction is " << v.transpose()
                << "\nThe surface normal is "
                << outward_surface_normal_vec.transpose();

            EXPECT_DOUBLE_EQ(0, bulk_normal.dot(outward_surface_normal_vec))
                << "The surface normal is not perpendicular to the element's "
                   "bulk normal."
                   "\nThe face id is "
                << face_id << "\nThe element's bulk normal is "
                << bulk_normal.transpose()
                << "\nThe surface normal of the edge is "
                << outward_surface_normal_vec.transpose();
        }
    }
    else if constexpr (dim == 3)
    {
        for (std::size_t face_id = 0; face_id < n_faces; ++face_id)
        {
            auto const face_ptr = ElementUtils::getFace(e, face_id);

            auto const outward_surface_normal_vec =
                calculateNormalizedSurfaceNormal(*face_ptr, e);

            for (std::size_t edge_id = 0;
                 edge_id < ElementUtils::getNumberOfFaces(*face_ptr);
                 ++edge_id)
            {
                auto const edge_ptr = ElementUtils::getFace(*face_ptr, edge_id);

                auto const& a = edge_ptr->getNode(0)->asEigenVector3d();
                auto const& b = edge_ptr->getNode(1)->asEigenVector3d();
                Eigen::Vector3d const v = b - a;

                EXPECT_DOUBLE_EQ(0, v.dot(outward_surface_normal_vec))
                    << "The normal is not perpendicular to the edge."
                       "\nThe face id is "
                    << face_id << "\nThe edge id is" << edge_id
                    << "\nThe edge's direction is " << v.transpose()
                    << "\nThe surface normal is "
                    << outward_surface_normal_vec.transpose();
            }
        }
    }
}

TYPED_TEST(MeshLibSurfaceNormalsTest, 2dReferenceElementNormalPointsToMinusZ)
{
    using MeshElementType = TypeParam;
    if constexpr (MeshElementType::dimension != 2)
    {
        return;
    }

    auto const& e = this->bulk_element;
    auto const n = MeshLib::FaceRule::getSurfaceNormal(e);

    // The fact that surface normals point in -z direction is the current OGS
    // convention. Usually one would expect the normals to point in +z
    // direction.
    EXPECT_DOUBLE_EQ(0, n[0])
        << "The surface normal does not point to -z but to " << n.transpose();
    EXPECT_DOUBLE_EQ(0, n[1])
        << "The surface normal does not point to -z but to " << n.transpose();
    EXPECT_GT(0, n[2]) << "The surface normal does not point to -z but to "
                       << n.transpose();
}
