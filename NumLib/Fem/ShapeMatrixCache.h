/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <boost/mp11.hpp>

#include "NumLib/Fem/FiniteElement/ElementTraitsLagrange.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"

namespace NumLib
{
namespace detail
{
template <typename ElementTraitsLagrange>
using GetMeshElement = typename ElementTraitsLagrange::Element;

template <typename ElementTraitsLagrange>
using GetShapeFunctionHigherOrder =
    typename ElementTraitsLagrange::ShapeFunction;

template <typename ElementTraitsLagrange>
using GetShapeFunctionLowerOrder =
    typename ElementTraitsLagrange::LowerOrderShapeFunction;

template <typename ShapeFunction>
using GetShapeMatrixPolicy =
    // dim does not matter, here, but
    // NumLib::detail::naturalCoordinatesMappingComputeShapeMatrices<>()
    // must be instantiated for this dim.
    ShapeMatrixPolicyType<ShapeFunction, 3 /* dim */>;

template <typename ShapeMatrixPolicy>
using GetShapeMatrix_N = typename ShapeMatrixPolicy::ShapeMatrices::ShapeType;
}  // namespace detail

class ShapeMatrixCache
{
    using ETLs = NumLib::AllElementTraitsLagrange;

    using MeshElements =
        boost::mp11::mp_transform<detail::GetMeshElement, ETLs>;

    // Higher order shape functions

    using ShapeFunctionsHigherOrder =
        boost::mp11::mp_transform<detail::GetShapeFunctionHigherOrder, ETLs>;
    using ShapeMatrixPoliciesHigherOrder =
        boost::mp11::mp_transform<detail::GetShapeMatrixPolicy,
                                  ShapeFunctionsHigherOrder>;
    using ShapeMatricesHigherOrder_N =
        boost::mp11::mp_transform<detail::GetShapeMatrix_N,
                                  ShapeMatrixPoliciesHigherOrder>;

    // std::tuple<std::vector<ShapeMatrix_N>, ...>
    using ShapeMatrixVectorsHigherOrder_N =
        boost::mp11::mp_transform<std::vector, ShapeMatricesHigherOrder_N>;

    static_assert(
        std::is_same_v<ShapeMatrixVectorsHigherOrder_N,
                       boost::mp11::mp_rename<ShapeMatrixVectorsHigherOrder_N,
                                              std::tuple>>,
        "The type alias ShapeMatrixVectorsHigherOrder_N must be a "
        "std::tuple<...>.");

    // Lower order shape functions

    using ShapeFunctionsLowerOrder =
        boost::mp11::mp_transform<detail::GetShapeFunctionLowerOrder, ETLs>;
    using ShapeMatrixPoliciesLowerOrder =
        boost::mp11::mp_transform<detail::GetShapeMatrixPolicy,
                                  ShapeFunctionsLowerOrder>;
    using ShapeMatricesLowerOrder_N =
        boost::mp11::mp_transform<detail::GetShapeMatrix_N,
                                  ShapeMatrixPoliciesLowerOrder>;

    // std::tuple<std::vector<ShapeMatrix_N>, ...>
    using ShapeMatrixVectorsLowerOrder_N =
        boost::mp11::mp_transform<std::vector, ShapeMatricesLowerOrder_N>;

    static_assert(
        std::is_same_v<
            ShapeMatrixVectorsLowerOrder_N,
            boost::mp11::mp_rename<ShapeMatrixVectorsLowerOrder_N, std::tuple>>,
        "The type alias ShapeMatrixVectorsLowerOrder_N must be a "
        "std::tuple<...>.");

public:
    explicit ShapeMatrixCache(unsigned const integration_order);

    template <typename MeshElement>
    auto const& NsHigherOrder() const
    {
        using Index = boost::mp11::mp_find<MeshElements, MeshElement>;
        return std::get<Index::value>(Nss_higher_order_);
    }

    template <typename MeshElement>
    auto const& NsLowerOrder() const
    {
        using Index = boost::mp11::mp_find<MeshElements, MeshElement>;
        return std::get<Index::value>(Nss_lower_order_);
    }

    /// Returns the number of elements in the ShapeMatrixVectorsHigherOrder_N
    /// tuple.
    static constexpr std::size_t size()
    {
        return boost::mp11::mp_size<ETLs>::value;
    }

    template <typename MeshElement>
    using ShapeFunctionHigherOrder =
        boost::mp11::mp_at<ShapeFunctionsHigherOrder,
                           boost::mp11::mp_find<MeshElements, MeshElement>>;

    template <typename MeshElement>
    using ShapeFunctionLowerOrder =
        boost::mp11::mp_at<ShapeFunctionsLowerOrder,
                           boost::mp11::mp_find<MeshElements, MeshElement>>;

private:
    ShapeMatrixVectorsHigherOrder_N Nss_higher_order_;
    ShapeMatrixVectorsLowerOrder_N Nss_lower_order_;
};
}  // namespace NumLib
