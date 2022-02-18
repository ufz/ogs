/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <gtest/gtest.h>

#include <unordered_map>

#include "MeshLib/Elements/Elements.h"
#include "NumLib/Fem/Integration/GaussLegendreIntegrationPolicy.h"
#include "NumLib/Fem/Integration/IntegrationMethodRegistry.h"

namespace
{
// helper type for converting lists of types
template <typename InputList, template <typename...> typename NewListType>
struct ConvertListType;
template <template <typename...> typename OldListType,
          typename... Ts,
          template <typename...>
          typename NewListType>
struct ConvertListType<OldListType<Ts...>, NewListType>
{
    using type = NewListType<Ts...>;
};
template <typename InputList, template <typename...> typename NewListType>
using ConvertListType_t =
    typename ConvertListType<InputList, NewListType>::type;

template <typename IntegrationMethod>
static unsigned getMaximumIntegrationOrderOfTemplatedIntegrationMethod(
    unsigned const order_cutoff)
{
    unsigned order = 1;
    int num_int_pts_last_iteration = -1;

    // trying all integration orders until one fails
    for (; order < order_cutoff; ++order)
    {
        int num_int_pts = -1;
        try
        {
            IntegrationMethod int_meth{order};
            num_int_pts = int_meth.getNumberOfPoints();
            unsigned const int_pt = 0;  // integration point zero should work
                                        // for any valid integration order
            int_meth.getWeightedPoint(int_pt);
        }
        catch (...)
        {
            // Integration method or weighted point is not available, so we have
            // exceeded the maximum integration order.
            --order;
            break;
        }

        if (num_int_pts <= num_int_pts_last_iteration)
        {
            // The number of integration points is not greater than for the
            // preceding integration order, so we have exceeded the maximum
            // integration order.
            --order;
            break;
        }

        num_int_pts_last_iteration = num_int_pts;
    }

    EXPECT_NE(order_cutoff, order)
        << "We hit the order cutoff. Determination of the maximum integration "
           "order in this unit test is wrong";

    return order;
}

template <typename MeshElementType>
static unsigned getMaximumIntegrationOrderFromIntegrationMethodRegistry(
    unsigned const order_cutoff)
{
    unsigned order = 1;

    // trying all integration orders until one fails
    for (; order < order_cutoff; ++order)
    {
        try
        {
            NumLib::IntegrationMethodRegistry::getIntegrationMethod<
                MeshElementType>(order);
        }
        catch (...)
        {
            // Integration method or weighted point is not available, so we have
            // exceeded the maximum integration order.
            --order;
            break;
        }
    }

    EXPECT_NE(order_cutoff, order)
        << "We hit the order cutoff. Determination of the maximum integration "
           "order in this unit test is wrong";

    return order;
}
}  // namespace

namespace MathLib
{
// Used to output helpful information if the unit tests fail.
static std::ostream& operator<<(std::ostream& os,
                                MathLib::WeightedPoint const& wp)
{
    auto const dim = wp.getDimension();
    os << "WP[" << dim << "D]{{";
    for (decltype(+dim) comp = 0; comp < 3; ++comp)
    {
        if (comp != 0)
        {
            os << ", ";
        }
        os << wp[comp];
    }
    os << "}, weight=" << wp.getWeight() << '}';
    return os;
}
}  // namespace MathLib

static std::unordered_map<std::type_index, unsigned> initMapTypeToMaxOrder()
{
    using namespace MeshLib;

    std::unordered_map<std::type_index, unsigned> map_type_to_max_order;

    auto add = [&](std::type_info const& info, unsigned order)
    {
        auto const [it, inserted] =
            map_type_to_max_order.emplace(std::type_index(info), order);
        if (!inserted)
        {
            throw std::runtime_error(
                "Maximum integration order already present for the "
                "current mesh element type. See file " __FILE__);
        }
    };

    add(typeid(Point), 4);
    add(typeid(Line), 4);
    add(typeid(Line3), 4);
    add(typeid(Quad), 4);
    add(typeid(Quad8), 4);
    add(typeid(Quad9), 4);
    add(typeid(Tri), 4);
    add(typeid(Tri6), 4);
    add(typeid(Hex), 4);
    add(typeid(Hex20), 4);
    add(typeid(Tet), 4);
    add(typeid(Tet10), 4);
    add(typeid(Prism), 4);
    add(typeid(Prism15), 4);
    add(typeid(Pyramid), 3);
    add(typeid(Pyramid13), 3);

    if (std::tuple_size_v<AllElementTypes> != map_type_to_max_order.size())
    {
        throw std::runtime_error(
            "Some element types are missing in the mapping from mesh element "
            "type to maximum integration order. See file " __FILE__);
    }

    return map_type_to_max_order;
}

static const std::unordered_map<std::type_index, unsigned>
    map_type_to_max_order = initMapTypeToMaxOrder();

// Test fixture.
template <typename MeshElementType>
class NumLibIntegrationMethodRegistryTest : public ::testing::Test
{
    std::type_index const type_index = std::type_index(typeid(MeshElementType));

protected:
    unsigned const max_order = map_type_to_max_order.at(type_index);

    NumLib::GenericIntegrationMethod const& getGenericIntegrationMethod(
        unsigned order) const
    {
        return NumLib::IntegrationMethodRegistry::getIntegrationMethod<
            MeshElementType>(order);
    }
};

using MeshElementTypes =
    ConvertListType_t<MeshLib::AllElementTypes, ::testing::Types>;

TYPED_TEST_SUITE(NumLibIntegrationMethodRegistryTest, MeshElementTypes);

TYPED_TEST(NumLibIntegrationMethodRegistryTest, HasRightOrder)
{
    for (unsigned order = 1; order <= this->max_order; ++order)
    {
        auto const& generic_int_meth = this->getGenericIntegrationMethod(order);

        EXPECT_EQ(order, generic_int_meth.getIntegrationOrder());
    }
}

TYPED_TEST(NumLibIntegrationMethodRegistryTest, HasRightIntegrationPoints)
{
    using MeshElementType = TypeParam;

    for (unsigned order = 1; order <= this->max_order; ++order)
    {
        auto const& generic_int_meth = this->getGenericIntegrationMethod(order);

        using IntegrationPolicy =
            NumLib::GaussLegendreIntegrationPolicy<MeshElementType>;
        using NonGenericIntegrationMethod =
            typename IntegrationPolicy::IntegrationMethod;

        NonGenericIntegrationMethod int_meth{order};

        auto const num_int_pts = generic_int_meth.getNumberOfPoints();
        auto const num_int_pts_expected = int_meth.getNumberOfPoints();
        ASSERT_EQ(num_int_pts_expected, num_int_pts)
            << "Wrong number of integration points for integration order "
            << order << '.';

        for (unsigned ip = 0; ip < num_int_pts; ++ip)
        {
            auto const& wp_expected = int_meth.getWeightedPoint(ip);
            auto const& wp_actual = generic_int_meth.getWeightedPoint(ip);
            EXPECT_EQ(wp_expected, wp_actual)
                << "Wrong integration point (#" << ip
                << ") for integration order " << order << '.'
                << "\nexpected: " << wp_expected << "\nactual:   " << wp_actual;
        }
    }
}

TYPED_TEST(NumLibIntegrationMethodRegistryTest, OrderZeroForbidden)
{
    unsigned const order = 0;

    ASSERT_ANY_THROW(this->getGenericIntegrationMethod(order));
}

// Assert that we don't forget any integration method in this test suite.
TYPED_TEST(NumLibIntegrationMethodRegistryTest, CheckWeTestUpToMaxOrder)
{
    using MeshElementType = TypeParam;

    if constexpr (std::is_same_v<MeshElementType, MeshLib::Point>)
    {
        // Points are special. For them all integration schemes for all
        // integration orders are the same.
        return;
    }

    using IntegrationPolicy =
        NumLib::GaussLegendreIntegrationPolicy<MeshElementType>;
    using NonGenericIntegrationMethod =
        typename IntegrationPolicy::IntegrationMethod;

    unsigned const cutoff = 100;  // Way beyond everything we ever expect.
    unsigned const order =
        getMaximumIntegrationOrderOfTemplatedIntegrationMethod<
            NonGenericIntegrationMethod>(cutoff);

    ASSERT_EQ(this->max_order, order)
        << "The maximum integration order used in this unit test suite differs "
           "from the (apparent) maximum integration order of the current "
           "integration method. Either we forgot to set the proper integration "
           "order in this unit test suite or the determination of the maximum "
           "integration order in the current unit test is wrong.";
}

// Assert that the integration method registry contains all integration orders
// of the underlying integration methods.
TYPED_TEST(NumLibIntegrationMethodRegistryTest, MaxOrderFromRegistry)
{
    using MeshElementType = TypeParam;

    if constexpr (std::is_same_v<MeshElementType, MeshLib::Prism> ||
                  std::is_same_v<MeshElementType, MeshLib::Prism15>)
    {
        GTEST_SKIP() << "Prisms currently behave rather special and cannot be "
                        "tested with the logic in this test.";
    }

    unsigned const cutoff = 100;  // Way beyond everything we ever expect.
    unsigned const order =
        getMaximumIntegrationOrderFromIntegrationMethodRegistry<
            MeshElementType>(cutoff);

    ASSERT_EQ(this->max_order, order)
        << "The maximum integration order tested in this unit test suite "
           "differs from the maximum integration order deposited in the "
           "integration method registry for the current mesh element type";
}
