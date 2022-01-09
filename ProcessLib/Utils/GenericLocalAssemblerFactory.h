/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
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
#include "NumLib/Fem/Integration/GaussLegendreIntegrationPolicy.h"

namespace ProcessLib
{
/// A functor creating a local assembler with shape functions corresponding to
/// the mesh element type.
template <typename LocalAssemblerInterface, typename... ConstructorArgs>
struct GenericLocalAssemblerFactory
{
    using LocAsmIntfPtr = std::unique_ptr<LocalAssemblerInterface>;
    using LocAsmBuilder =
        std::function<LocAsmIntfPtr(MeshLib::Element const& e,
                                    std::size_t const local_matrix_size,
                                    ConstructorArgs&&...)>;

protected:  // only allow instances of subclasses
    explicit GenericLocalAssemblerFactory(
        NumLib::LocalToGlobalIndexMap const& dof_table)
        : _dof_table(dof_table)
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

protected:
    /// Mapping of element types to local assembler builders.
    std::unordered_map<std::type_index, LocAsmBuilder> _builders;
};

template <typename ShapeFunction, typename LocalAssemblerInterface,
          template <typename /* shp fct */, typename /* int meth */,
                    int /* global dim */>
          class LocalAssemblerImplementation,
          int GlobalDim, typename... ConstructorArgs>
class LocalAssemblerBuilderFactory
{
    using GLAF = GenericLocalAssemblerFactory<LocalAssemblerInterface,
                                              ConstructorArgs...>;
    using LocAsmIntfPtr = typename GLAF::LocAsmIntfPtr;
    using LocAsmBuilder = typename GLAF::LocAsmBuilder;

    using IntegrationMethod = typename NumLib::GaussLegendreIntegrationPolicy<
        typename ShapeFunction::MeshElement>::IntegrationMethod;

    using LocAsmImpl =
        LocalAssemblerImplementation<ShapeFunction, IntegrationMethod,
                                     GlobalDim>;

    LocalAssemblerBuilderFactory() = delete;

public:
    /// Generates a function that creates a new local assembler of type
    /// \c LocAsmImpl.
    static LocAsmBuilder create()
    {
        return [](MeshLib::Element const& e,
                  std::size_t const local_matrix_size,
                  ConstructorArgs&&... args)
        {
            return std::make_unique<LocAsmImpl>(
                e, local_matrix_size, std::forward<ConstructorArgs>(args)...);
        };
    }
};

}  // namespace ProcessLib
