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

#include <functional>
#include <memory>
#include <typeindex>
#include <unordered_map>

#include "NumLib/DOF/LocalToGlobalIndexMap.h"
#include "NumLib/Fem/Integration/IntegrationMethodProvider.h"

namespace ProcessLib
{
/// Used for template type checks with createLocalAssemblers... functions.
///
/// With this concept the createLocalAssemblers... functions transparently
/// support both OGS's default integration methods (requested via a
/// NumLib::IntegrationOrder argument) and arbitrary
/// NumLib::IntegrationMethodProvider's.
template <typename T>
concept IntegrationMethodProviderOrIntegrationOrder =
    NumLib::IntegrationMethodProvider<T> ||
    std::same_as<T, NumLib::IntegrationOrder>;

/// A functor creating a local assembler with shape functions corresponding to
/// the mesh element type.
template <typename LocalAssemblerInterface,
          NumLib::IntegrationMethodProvider IntegrationMethodProvider,
          typename... ConstructorArgs>
struct GenericLocalAssemblerFactory
{
    using LocAsmIntfPtr = std::unique_ptr<LocalAssemblerInterface>;
    using LocAsmBuilder = std::function<LocAsmIntfPtr(
        MeshLib::Element const& e,
        std::size_t const local_matrix_size,
        IntegrationMethodProvider const& integration_method_provider,
        ConstructorArgs&&...)>;

protected:  // only allow instances of subclasses
    explicit GenericLocalAssemblerFactory(
        NumLib::LocalToGlobalIndexMap const& dof_table,
        IntegrationMethodProvider const& integration_method_provider)
        : _dof_table{dof_table},
          _integration_method_provider{integration_method_provider}
    {
    }

public:
    /// Returns the newly created local assembler.
    ///
    /// \attention
    /// The index \c id is not necessarily the mesh item's id. Especially when
    /// having multiple meshes it will differ from the latter.
    LocAsmIntfPtr operator()(std::size_t const id,
                             MeshLib::Element const& mesh_item,
                             ConstructorArgs&&... args) const
    {
        auto const type_idx = std::type_index(typeid(mesh_item));
        auto const it = _builders.find(type_idx);

        if (it != _builders.end())
        {
            auto const num_local_dof = _dof_table.getNumberOfElementDOF(id);
            return it->second(mesh_item, num_local_dof,
                              _integration_method_provider,
                              std::forward<ConstructorArgs>(args)...);
        }
        OGS_FATAL(
            "You are trying to build a local assembler for an unknown mesh "
            "element type ({:s})."
            " Maybe you have disabled this mesh element type in your build "
            "configuration, or a mesh element order does not match shape "
            "function order given in the project file.",
            type_idx.name());
    }

private:
    NumLib::LocalToGlobalIndexMap const& _dof_table;
    IntegrationMethodProvider const& _integration_method_provider;

protected:
    /// Mapping of element types to local assembler builders.
    std::unordered_map<std::type_index, LocAsmBuilder> _builders;
};

template <typename ShapeFunction, typename LocalAssemblerInterface,
          template <typename /* shp fct */, int /* global dim */>
          class LocalAssemblerImplementation,
          NumLib::IntegrationMethodProvider IntegrationMethodProvider,
          int GlobalDim, typename... ConstructorArgs>
class LocalAssemblerBuilderFactory
{
    using GLAF = GenericLocalAssemblerFactory<
        LocalAssemblerInterface, IntegrationMethodProvider, ConstructorArgs...>;
    using LocAsmIntfPtr = typename GLAF::LocAsmIntfPtr;
    using LocAsmBuilder = typename GLAF::LocAsmBuilder;

    using LocAsmImpl = LocalAssemblerImplementation<ShapeFunction, GlobalDim>;

    LocalAssemblerBuilderFactory() = delete;

public:
    /// Generates a function that creates a new local assembler of type
    /// \c LocAsmImpl.
    template <typename MeshElement>
    static LocAsmBuilder create()
    {
        return [](MeshLib::Element const& e,
                  std::size_t const local_matrix_size,
                  IntegrationMethodProvider const& integration_method_provider,
                  ConstructorArgs&&... args)
        {
            auto const& integration_method =
                integration_method_provider
                    .template getIntegrationMethod<MeshElement>(e);

            static_assert(
                std::is_constructible_v<
                    LocAsmImpl, MeshLib::Element const&, std::size_t,
                    decltype(integration_method), ConstructorArgs&&...>,
                "The given local assembler implementation is not "
                "constructible from the provided arguments.");

            return std::make_unique<LocAsmImpl>(
                e, local_matrix_size, integration_method,
                std::forward<ConstructorArgs>(args)...);
        };
    }
};

}  // namespace ProcessLib
