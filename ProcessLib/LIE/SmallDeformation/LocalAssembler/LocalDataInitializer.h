/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <functional>
#include <memory>
#include <type_traits>
#include <typeindex>
#include <typeinfo>
#include <unordered_map>

#include "MeshLib/Elements/Elements.h"
#include "NumLib/DOF/LocalToGlobalIndexMap.h"
#include "NumLib/Fem/Integration/GaussLegendreIntegrationPolicy.h"
#include "ProcessLib/Utils/EnabledElements.h"

namespace ProcessLib
{
namespace LIE
{
namespace SmallDeformation
{
/// The LocalDataInitializer is a functor creating a local assembler data with
/// corresponding to the mesh element type shape functions and calling
/// initialization of the new local assembler data.
/// For example for MeshLib::Quad a local assembler data with template argument
/// NumLib::ShapeQuad4 is created.
template <typename LocalAssemblerInterface,
          template <typename, typename, int> class LocalAssemblerDataMatrix,
          template <typename, typename, int>
          class LocalAssemblerDataMatrixNearFracture,
          template <typename, typename, int> class LocalAssemblerDataFracture,
          int GlobalDim, typename... ConstructorArgs>
class LocalDataInitializer final
{
    struct IsElementEnabled
    {
        template <typename ElementTraits>
        constexpr bool operator()(ElementTraits*) const
        {
            if constexpr (GlobalDim < ElementTraits::ShapeFunction::DIM)
            {
                return false;
            }

            // exclude 0D elements
            return ElementTraits::Element::dimension >= 1;
        }
    };

public:
    using LADataIntfPtr = std::unique_ptr<LocalAssemblerInterface>;

    explicit LocalDataInitializer(
        NumLib::LocalToGlobalIndexMap const& dof_table)
        : _dof_table(dof_table)
    {
        // REMARKS: At the moment, only a 2D mesh with 1D elements are
        // supported.

        using EnabledElementTraits =
            decltype(BaseLib::TMP::filter<EnabledElementTraitsLagrange>(
                std::declval<IsElementEnabled>()));

        BaseLib::TMP::foreach<EnabledElementTraits>(
            [this]<typename ET>(ET*)
            {
                using Elt = typename ET::Element;
                using Shp = typename ET::ShapeFunction;

                _builder[std::type_index(typeid(Elt))] =
                    makeLocalAssemblerBuilder<Shp>();
            });
    }

    /// Returns data pointer to the newly created local assembler data.
    ///
    /// \attention
    /// The index \c id is not necessarily the mesh item's id. Especially when
    /// having multiple meshes it will differ from the latter.
    LADataIntfPtr operator()(std::size_t const id,
                    MeshLib::Element const& mesh_item,
                    ConstructorArgs&&... args) const
    {
        auto const type_idx = std::type_index(typeid(mesh_item));
        auto const it = _builder.find(type_idx);

        if (it == _builder.end())
        {
            OGS_FATAL(
                "You are trying to build a local assembler for an unknown mesh "
                "element type ({:s})."
                " Maybe you have disabled this mesh element type in your build "
                "configuration, or a mesh element order does not match shape "
                "function order given in the project file.",
                type_idx.name());
        }

        auto const n_local_dof = _dof_table.getNumberOfElementDOF(id);
        auto const n_global_components =
            _dof_table.getNumberOfElementComponents(id);
        auto const varIDs = _dof_table.getElementVariableIDs(id);

        std::vector<unsigned> dofIndex_to_localIndex;
        if (mesh_item.getDimension() < GlobalDim ||
            n_global_components > GlobalDim)
        {
            dofIndex_to_localIndex.resize(n_local_dof);
            unsigned dof_id = 0;
            unsigned local_id = 0;
            for (auto i : varIDs)
            {
                for (int j = 0; j < _dof_table.getNumberOfVariableComponents(i);
                     j++)
                {
                    auto const& ms = _dof_table.getMeshSubset(i, j);
                    auto const mesh_id = ms.getMeshID();
                    for (unsigned k = 0; k < mesh_item.getNumberOfNodes(); k++)
                    {
                        MeshLib::Location l(mesh_id,
                                            MeshLib::MeshItemType::Node,
                                            getNodeIndex(mesh_item, k));
                        auto global_index = _dof_table.getGlobalIndex(l, i, j);
                        if (global_index != NumLib::MeshComponentMap::nop)
                        {
                            dofIndex_to_localIndex[dof_id++] = local_id;
                        }
                        local_id++;
                    }
                }
            }
        }

        return it->second(mesh_item, varIDs.size(), n_local_dof,
                          dofIndex_to_localIndex,
                          std::forward<ConstructorArgs>(args)...);
    }

private:
    using LADataBuilder = std::function<LADataIntfPtr(
        MeshLib::Element const& e,
        std::size_t const n_variables,
        std::size_t const local_matrix_size,
        std::vector<unsigned> const& dofIndex_to_localIndex,
        ConstructorArgs&&...)>;

    template <typename ShapeFunction>
    using IntegrationMethod = typename NumLib::GaussLegendreIntegrationPolicy<
        typename ShapeFunction::MeshElement>::IntegrationMethod;

    template <typename ShapeFunction>
    using LADataMatrix =
        LocalAssemblerDataMatrix<ShapeFunction,
                                 IntegrationMethod<ShapeFunction>, GlobalDim>;

    /// Generates a function that creates a new LocalAssembler of type
    /// LAData<ShapeFunction>. Only functions with shape function's dimension
    /// less or equal to the global dimension are instantiated, e.g.  following
    /// combinations of shape functions and global dimensions: (Line2, 1),
    /// (Line2, 2), (Line2, 3), (Hex20, 3) but not (Hex20, 2) or (Hex20, 1).
    template <typename ShapeFunction>
    static LADataBuilder makeLocalAssemblerBuilder()
    {
        return [](MeshLib::Element const& e,
                  std::size_t const n_variables,
                  std::size_t const local_matrix_size,
                  std::vector<unsigned> const& dofIndex_to_localIndex,
                  ConstructorArgs&&... args)
        {
            if (e.getDimension() == GlobalDim)
            {
                if (dofIndex_to_localIndex.empty())
                {
                    return LADataIntfPtr{new LADataMatrix<ShapeFunction>{
                        e, local_matrix_size,
                        std::forward<ConstructorArgs>(args)...}};
                }

                return LADataIntfPtr{
                    new LADataMatrixNearFracture<ShapeFunction>{
                        e, n_variables, local_matrix_size,
                        dofIndex_to_localIndex,
                        std::forward<ConstructorArgs>(args)...}};
            }
            return LADataIntfPtr{new LAFractureData<ShapeFunction>{
                e, n_variables, local_matrix_size, dofIndex_to_localIndex,
                std::forward<ConstructorArgs>(args)...}};
        };
    }

    /// Mapping of element types to local assembler constructors.
    std::unordered_map<std::type_index, LADataBuilder> _builder;

    NumLib::LocalToGlobalIndexMap const& _dof_table;

    template <typename ShapeFunction>
    using LADataMatrixNearFracture = LocalAssemblerDataMatrixNearFracture<
        ShapeFunction, IntegrationMethod<ShapeFunction>, GlobalDim>;

    template <typename ShapeFunction>
    using LAFractureData =
        LocalAssemblerDataFracture<ShapeFunction,
                                   IntegrationMethod<ShapeFunction>, GlobalDim>;
};

}  // namespace SmallDeformation
}  // namespace LIE
}  // namespace ProcessLib
