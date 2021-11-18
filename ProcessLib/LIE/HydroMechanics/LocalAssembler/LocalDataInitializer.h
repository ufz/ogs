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
#include "NumLib/Fem/FiniteElement/LowerDimShapeTable.h"
#include "NumLib/Fem/Integration/GaussLegendreIntegrationPolicy.h"
#include "ProcessLib/Utils/EnabledElements.h"

namespace ProcessLib
{
namespace LIE
{
namespace HydroMechanics
{
/// The LocalDataInitializer is a functor creating a local assembler data with
/// corresponding to the mesh element type shape functions and calling
/// initialization of the new local assembler data.
/// For example for MeshLib::Quad a local assembler data with template argument
/// NumLib::ShapeQuad4 is created.
///
/// \attention This is modified version of the ProcessLib::LocalDataInitializer
/// class which does not include line elements, allows only shapefunction of
/// order 2, and provides additional template argument GlobalDim.
template <typename LocalAssemblerInterface,
          template <typename, typename, typename, int>
          class LocalAssemblerDataMatrix,
          template <typename, typename, typename, int>
          class LocalAssemblerDataMatrixNearFracture,
          template <typename, typename, typename, int>
          class LocalAssemblerDataFracture,
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

            // exclude 0D elements and linear elements
            return ElementTraits::Element::dimension >= 1 &&
                   ElementTraits::ShapeFunction::ORDER == 2;
        }
    };

public:
    using LADataIntfPtr = std::unique_ptr<LocalAssemblerInterface>;

    explicit LocalDataInitializer(
        NumLib::LocalToGlobalIndexMap const& dof_table)
        : _dof_table(dof_table)
    {
        using EnabledElementTraits =
            decltype(BaseLib::TMP::filter<EnabledElementTraitsLagrange>(
                std::declval<IsElementEnabled>()));

        BaseLib::TMP::foreach<EnabledElementTraits>(
            [this]<typename ET>(ET*)
            {
                using Elt = typename ET::Element;
                using Shp = typename ET::ShapeFunction;
                using ShpLow = typename ET::LowerOrderShapeFunction;

                _builder[std::type_index(typeid(Elt))] =
                    makeLocalAssemblerBuilder<Shp, ShpLow>();
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
        auto const varIDs = _dof_table.getElementVariableIDs(id);
        bool const isPressureDeactivated = (varIDs.front() != 0);
        std::vector<int> involved_varIDs;  // including deactivated elements
        involved_varIDs.reserve(varIDs.size() + 1);
        if (isPressureDeactivated)
        {
            involved_varIDs.push_back(0);  // always pressure come in
        }
        involved_varIDs.insert(involved_varIDs.end(), varIDs.begin(),
                               varIDs.end());

        std::vector<unsigned> dofIndex_to_localIndex;

        // matrix and fracture assemblers with enrichments
        dofIndex_to_localIndex.resize(n_local_dof);
        std::vector<unsigned> vec_n_element_nodes;
        // TODO how to get the shape function order for each variable?
        vec_n_element_nodes.push_back(
            mesh_item.getNumberOfBaseNodes());  // pressure
        auto const max_varID = *std::max_element(varIDs.begin(), varIDs.end());
        for (int i = 1; i < max_varID + 1; i++)
        {
            vec_n_element_nodes.push_back(
                mesh_item.getNumberOfNodes());  // displacements
        }

        unsigned local_id = 0;
        unsigned dof_id = 0;
        for (unsigned i = 0; i < involved_varIDs.size(); i++)
        {
            auto const var_id = involved_varIDs[i];
            auto const n_var_comp =
                _dof_table.getNumberOfVariableComponents(var_id);
            auto const n_var_element_nodes = vec_n_element_nodes[i];
            for (int var_comp_id = 0; var_comp_id < n_var_comp; var_comp_id++)
            {
                auto const& ms = _dof_table.getMeshSubset(var_id, var_comp_id);
                auto const mesh_id = ms.getMeshID();
                for (unsigned k = 0; k < n_var_element_nodes; k++)
                {
                    MeshLib::Location l(mesh_id,
                                        MeshLib::MeshItemType::Node,
                                        getNodeIndex(mesh_item, k));
                    auto global_index =
                        _dof_table.getGlobalIndex(l, var_id, var_comp_id);
                    if (global_index != NumLib::MeshComponentMap::nop)
                    {
                        dofIndex_to_localIndex[dof_id++] = local_id;
                    }
                    local_id++;
                }
            }
        }

        return it->second(mesh_item, involved_varIDs.size(), n_local_dof,
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

    template <typename ShapeFunctionDisplacement>
    using IntegrationMethod = typename NumLib::GaussLegendreIntegrationPolicy<
        typename ShapeFunctionDisplacement::MeshElement>::IntegrationMethod;

    template <typename ShapeFunctionDisplacement,
              typename ShapeFunctionPressure>
    using LADataMatrix = LocalAssemblerDataMatrix<
        ShapeFunctionDisplacement, ShapeFunctionPressure,
        IntegrationMethod<ShapeFunctionDisplacement>, GlobalDim>;

    template <typename ShapeFunctionDisplacement,
              typename ShapeFunctionPressure>
    using LADataMatrixNearFracture = LocalAssemblerDataMatrixNearFracture<
        ShapeFunctionDisplacement, ShapeFunctionPressure,
        IntegrationMethod<ShapeFunctionDisplacement>, GlobalDim>;

    template <typename ShapeFunctionDisplacement,
              typename ShapeFunctionPressure>
    using LAFractureData = LocalAssemblerDataFracture<
        ShapeFunctionDisplacement, ShapeFunctionPressure,
        IntegrationMethod<ShapeFunctionDisplacement>, GlobalDim>;

    /// Generates a function that creates a new LocalAssembler of type
    /// LAData<ShapeFunctionDisplacement>. Only functions with shape function's
    /// dimension less or equal to the global dimension are instantiated, e.g.
    /// following combinations of shape functions and global dimensions: (Line2,
    /// 1),
    /// (Line2, 2), (Line2, 3), (Hex20, 3) but not (Hex20, 2) or (Hex20, 1).
    template <typename ShapeFunctionDisplacement,
              typename ShapeFunctionPressure>
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
                if (n_variables == 2)
                {
                    return LADataIntfPtr{
                        new LADataMatrix<ShapeFunctionDisplacement,
                                         ShapeFunctionPressure>{
                            e, n_variables, local_matrix_size,
                            dofIndex_to_localIndex,
                            std::forward<ConstructorArgs>(args)...}};
                }
                return LADataIntfPtr{new LADataMatrixNearFracture<
                    ShapeFunctionDisplacement, ShapeFunctionPressure>{
                    e, n_variables, local_matrix_size, dofIndex_to_localIndex,
                    std::forward<ConstructorArgs>(args)...}};
            }
            return LADataIntfPtr{new LAFractureData<ShapeFunctionDisplacement,
                                                    ShapeFunctionPressure>{
                e, local_matrix_size, dofIndex_to_localIndex,
                std::forward<ConstructorArgs>(args)...}};
        };
    }

    /// Mapping of element types to local assembler constructors.
    std::unordered_map<std::type_index, LADataBuilder> _builder;

    NumLib::LocalToGlobalIndexMap const& _dof_table;
};

}  // namespace HydroMechanics
}  // namespace LIE
}  // namespace ProcessLib
