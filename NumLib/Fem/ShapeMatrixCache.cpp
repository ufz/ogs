/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "ShapeMatrixCache.h"

#include <range/v3/range/conversion.hpp>
#include <range/v3/view/transform.hpp>

#include "NumLib/Fem/InitShapeMatrices.h"
#include "NumLib/Fem/Integration/IntegrationMethodRegistry.h"
#include "NumLib/Fem/ReferenceElement.h"
namespace
{
template <typename ShapeFunction>
static auto initShapeMatrices(unsigned const integration_order,
                              boost::mp11::mp_identity<ShapeFunction>)
{
    using namespace boost::mp11;

    using MeshElement = typename ShapeFunction::MeshElement;

    auto const& integration_method =
        NumLib::IntegrationMethodRegistry::getIntegrationMethod<MeshElement>(
            NumLib::IntegrationOrder{integration_order});

    using ShapeMatrixPolicy =
        NumLib::detail::GetShapeMatrixPolicy<ShapeFunction>;

    // The reference element will be used as dummy argument below. That's
    // sufficient, since we only want to compute the shape matrix N.
    NumLib::ReferenceElement<MeshElement> reference_element;

    // dim does not matter, here.
    constexpr int dim = 3;
    auto sms = NumLib::initShapeMatrices<ShapeFunction,
                                         ShapeMatrixPolicy,
                                         dim,
                                         NumLib::ShapeMatrixType::N>(
        reference_element.element,
        false /*is_axially_symmetric*/,
        integration_method);

    return sms |
           ranges::views::transform([](auto& sm) { return std::move(sm.N); }) |
           ranges::to<std::vector>;
}
}  // namespace

namespace NumLib
{
ShapeMatrixCache::ShapeMatrixCache(unsigned const integration_order)
{
    using namespace boost::mp11;
    using Indices = mp_iota<mp_size<ETLs>>;

    mp_for_each<Indices>(
        [this, integration_order](auto index)
        {
            using ETL = mp_at_c<ETLs, index>;
            using ShapeFunctionHigherOrder =
                detail::GetShapeFunctionHigherOrder<ETL>;

            std::get<index.value>(Nss_higher_order_) = ::initShapeMatrices(
                integration_order, mp_identity<ShapeFunctionHigherOrder>{});

            using ShapeFunctionLowerOrder =
                detail::GetShapeFunctionLowerOrder<ETL>;

            std::get<index.value>(Nss_lower_order_) = ::initShapeMatrices(
                integration_order, mp_identity<ShapeFunctionLowerOrder>{});
        });
}

}  // namespace NumLib
