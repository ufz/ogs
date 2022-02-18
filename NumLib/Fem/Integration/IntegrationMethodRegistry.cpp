/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "IntegrationMethodRegistry.h"

#include <unordered_map>

#include "BaseLib/Error.h"
#include "GaussLegendreIntegrationPolicy.h"
#include "MeshLib/Elements/Elements.h"

static constexpr unsigned MAX_ORDER_REGULAR = 4u;

template <typename MeshElement>
using IntegrationPolicy = NumLib::GaussLegendreIntegrationPolicy<MeshElement>;

template <typename IntegrationPolicy_>
static NumLib::GenericIntegrationMethod createGenericIntegrationMethod(
    unsigned const order)
{
    typename IntegrationPolicy_::IntegrationMethod meth{order};
    unsigned const np = meth.getNumberOfPoints();

    std::vector<MathLib::WeightedPoint> wps;
    wps.reserve(np);

    for (unsigned ip = 0; ip < np; ++ip)
    {
        wps.emplace_back(meth.getWeightedPoint(ip));

        if (wps.back() != meth.getWeightedPoint(ip))
        {
            throw std::runtime_error(
                "createGenericIntegrationMethod mismatch for ip=" +
                std::to_string(ip) + ", order=" + std::to_string(order) +
                ", method=" + typeid(decltype(meth)).name());
        }
    }

    return NumLib::GenericIntegrationMethod{order, std::move(wps)};
}

template <typename MeshElement>
static void putIntegrationMethodsFor(
    std::unordered_map<std::type_index,
                       std::vector<NumLib::GenericIntegrationMethod>>&
        integration_methods_by_mesh_element_type,
    unsigned max_order)
{
    std::vector<NumLib::GenericIntegrationMethod> integration_methods;
    integration_methods.reserve(max_order + 1);

    // order 0 -> no valid weighted points
    NumLib::GenericIntegrationMethod invalidMethod{0, {}};
    integration_methods.push_back(std::move(invalidMethod));

    for (unsigned order = 1; order <= max_order; ++order)
    {
        integration_methods.push_back(
            createGenericIntegrationMethod<IntegrationPolicy<MeshElement>>(
                order));
    }

    integration_methods_by_mesh_element_type.emplace(
        std::type_index(typeid(MeshElement)), std::move(integration_methods));
}

static void putIntegrationMethodsForDim0(
    std::unordered_map<std::type_index,
                       std::vector<NumLib::GenericIntegrationMethod>>&
        integration_methods_by_mesh_element_type)
{
    putIntegrationMethodsFor<MeshLib::Point>(
        integration_methods_by_mesh_element_type,
        MAX_ORDER_REGULAR /* arbitrary cutoff */);
}

static void putIntegrationMethodsForDim1(
    std::unordered_map<std::type_index,
                       std::vector<NumLib::GenericIntegrationMethod>>&
        integration_methods_by_mesh_element_type)
{
    putIntegrationMethodsFor<MeshLib::Line>(
        integration_methods_by_mesh_element_type, MAX_ORDER_REGULAR);

    putIntegrationMethodsFor<MeshLib::Line3>(
        integration_methods_by_mesh_element_type, MAX_ORDER_REGULAR);
}

static void putIntegrationMethodsForDim2(
    std::unordered_map<std::type_index,
                       std::vector<NumLib::GenericIntegrationMethod>>&
        integration_methods_by_mesh_element_type)
{
    putIntegrationMethodsFor<MeshLib::Quad>(
        integration_methods_by_mesh_element_type, MAX_ORDER_REGULAR);

    putIntegrationMethodsFor<MeshLib::Quad8>(
        integration_methods_by_mesh_element_type, MAX_ORDER_REGULAR);

    putIntegrationMethodsFor<MeshLib::Quad9>(
        integration_methods_by_mesh_element_type, MAX_ORDER_REGULAR);

    putIntegrationMethodsFor<MeshLib::Tri>(
        integration_methods_by_mesh_element_type, 4);

    putIntegrationMethodsFor<MeshLib::Tri6>(
        integration_methods_by_mesh_element_type, 4);
}

static void putIntegrationMethodsForDim3(
    std::unordered_map<std::type_index,
                       std::vector<NumLib::GenericIntegrationMethod>>&
        integration_methods_by_mesh_element_type)
{
    putIntegrationMethodsFor<MeshLib::Hex>(
        integration_methods_by_mesh_element_type, MAX_ORDER_REGULAR);

    putIntegrationMethodsFor<MeshLib::Hex20>(
        integration_methods_by_mesh_element_type, MAX_ORDER_REGULAR);

    putIntegrationMethodsFor<MeshLib::Tet>(
        integration_methods_by_mesh_element_type, 4);

    putIntegrationMethodsFor<MeshLib::Tet10>(
        integration_methods_by_mesh_element_type, 4);

    // TODO max order 2? cf. IntegrationGaussLegendrePrism::getNumberOfPoints
    // (only order 2)
    putIntegrationMethodsFor<MeshLib::Prism>(
        integration_methods_by_mesh_element_type, 2);

    putIntegrationMethodsFor<MeshLib::Prism15>(
        integration_methods_by_mesh_element_type, 2);

    putIntegrationMethodsFor<MeshLib::Pyramid>(
        integration_methods_by_mesh_element_type, 3);

    putIntegrationMethodsFor<MeshLib::Pyramid13>(
        integration_methods_by_mesh_element_type, 3);
}

static std::unordered_map<std::type_index,
                          std::vector<NumLib::GenericIntegrationMethod>>
initIntegrationMethods()
{
    std::unordered_map<std::type_index,
                       std::vector<NumLib::GenericIntegrationMethod>>
        integration_methods_by_mesh_element_type;

    putIntegrationMethodsForDim0(integration_methods_by_mesh_element_type);
    putIntegrationMethodsForDim1(integration_methods_by_mesh_element_type);
    putIntegrationMethodsForDim2(integration_methods_by_mesh_element_type);
    putIntegrationMethodsForDim3(integration_methods_by_mesh_element_type);

    return integration_methods_by_mesh_element_type;
}

namespace NumLib::IntegrationMethodRegistry
{
GenericIntegrationMethod const& getIntegrationMethod(
    std::type_index const mesh_element_type, unsigned const order)
{
    if (order == 0)
    {
        // For 0D elements (points) the order does not matter, still we don't
        // want order zero there. For all other element types integration order
        // does matter.
        OGS_FATAL("An integration order of 0 is not supported.");
    }

    // The actual integration method registry.
    //
    // A mapping \code [mesh element type] -> [list of integration
    // methods]\endcode, where each entry in the list corresponds to integration
    // order 0, 1, 2, ...
    //
    // Having this as a static __local__ variable circumvents the static
    // initialization order fiasco we would have if this was a static __global__
    // variable. The fiasco arises, because this registry uses static global
    // data from MathLib. Currently (Feb 2022) we cannot circumvent the problem
    // with constinit due to lack of compiler support.
    static const std::unordered_map<
        std::type_index,
        std::vector<NumLib::GenericIntegrationMethod>>
        integration_methods_by_mesh_element_type = initIntegrationMethods();

    if (auto it =
            integration_methods_by_mesh_element_type.find(mesh_element_type);
        it != integration_methods_by_mesh_element_type.end())
    {
        auto& integration_methods = it->second;

        if (order >= integration_methods.size())
        {
            OGS_FATAL(
                "Integration order {} is not supported for mesh elements of "
                "type {}",
                order,
                mesh_element_type.name());
        }

        return integration_methods[order];
    }

    OGS_FATAL(
        "No integration methods are available for mesh elements of type {}",
        mesh_element_type.name());
}
}  // namespace NumLib::IntegrationMethodRegistry
