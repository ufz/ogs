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

#include "IntegrationMethodRegistry.h"
#include "MeshLib/Elements/Quad.h"

namespace NumLib
{
/// An \c IntegrationMethodProvider provides <tt>GenericIntegrationMethod</tt>s
/// for <tt>MeshLib::Element</tt>s.
template <typename T>
concept IntegrationMethodProvider = requires(T t)
{
    {
        // using Quad as a representative mesh element type
        t.template getIntegrationMethod<MeshLib::Quad>(
            std::declval<MeshLib::Element>())
        } -> std::same_as<GenericIntegrationMethod const&>;
};

/// Provides OGS's default Gauss-Legendre integration method of a given order
/// for any mesh element.
class DefaultIntegrationMethodProvider
{
public:
    explicit DefaultIntegrationMethodProvider(
        IntegrationOrder const integration_order)
        : integration_order_{integration_order}
    {
    }

    template <typename MeshElement>
    GenericIntegrationMethod const& getIntegrationMethod(
        MeshLib::Element const& /*e*/) const
    {
        return IntegrationMethodRegistry::template getIntegrationMethod<
            MeshElement>(integration_order_);
    }

private:
    IntegrationOrder integration_order_;
};

/// Used in template functions to get the \c DefaultIntegrationMethodProvider
/// for the given \c integration_order.
inline DefaultIntegrationMethodProvider getIntegrationMethodProvider(
    NumLib::IntegrationOrder const integration_order)
{
    return DefaultIntegrationMethodProvider{integration_order};
}

}  // namespace NumLib
