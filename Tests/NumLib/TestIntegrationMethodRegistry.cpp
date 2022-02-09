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
    add(typeid(Prism), 1);  // Actually, there is only one unique integration
                            // method for prisms, currently.
    add(typeid(Prism15), 1);
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
