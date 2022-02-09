/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <typeindex>

#include "GenericIntegrationMethod.h"

namespace NumLib
{
namespace IntegrationMethodRegistry
{
//! Provides the default integration method for a given mesh element type and a
//! given integration order.
GenericIntegrationMethod const& getIntegrationMethod(
    std::type_index const mesh_element_type, unsigned const order);

//! Templated version of above method. Provided for convenience.
template <typename MeshElement>
GenericIntegrationMethod const& getIntegrationMethod(unsigned order)
{
    return getIntegrationMethod(std::type_index(typeid(MeshElement)), order);
}
}  // namespace IntegrationMethodRegistry
}  // namespace NumLib
