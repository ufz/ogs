// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <typeindex>

#include "GenericIntegrationMethod.h"

namespace NumLib
{
//! A named type avoiding that the integration order is confused with other
//! integral types when passed as function argument.
struct IntegrationOrder
{
    explicit IntegrationOrder(unsigned const order_) : order{order_} {}

    unsigned order;
};

namespace IntegrationMethodRegistry
{
//! Provides the default integration method for a given mesh element type and a
//! given integration order.
GenericIntegrationMethod const& getIntegrationMethod(
    std::type_index const mesh_element_type, IntegrationOrder const order);

//! Templated version of above method. Provided for convenience.
template <typename MeshElement>
GenericIntegrationMethod const& getIntegrationMethod(
    IntegrationOrder const order)
{
    return getIntegrationMethod(std::type_index(typeid(MeshElement)), order);
}
}  // namespace IntegrationMethodRegistry
}  // namespace NumLib
