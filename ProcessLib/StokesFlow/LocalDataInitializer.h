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

namespace ProcessLib::StokesFlow
{
/// The LocalDataInitializer is a functor creating a local assembler data with
/// corresponding to the mesh element type shape functions and calling
/// initialization of the new local assembler data.
/// For example for MeshLib::Quad a local assembler data with template argument
/// NumLib::ShapeQuad4 is created.
///
/// \attention This is modified version of the ProcessLib::LocalDataInitializer
/// class which does not include line or point elements. For the shape functions
/// of order 2 (used for fluid velocity in StokesFlow) a shape function of order
/// 1 will be used for the scalar variables.
template <typename LocalAssemblerInterface,
          template <typename, typename, typename, int> class LocalAssemblerData,
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

            return ElementTraits::Element::dimension >= 2 &&
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
                using MeshElement = typename ET::Element;
                using ShapeFunction = typename ET::ShapeFunction;
                using LowerOrderShapeFunction =
                    typename ET::LowerOrderShapeFunction;

                _builder[std::type_index(typeid(MeshElement))] =
                    makeLocalAssemblerBuilder<ShapeFunction,
                                              LowerOrderShapeFunction>();
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

        if (it != _builder.end())
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
    using LADataBuilder =
        std::function<LADataIntfPtr(MeshLib::Element const& e,
                                    std::size_t const local_matrix_size,
                                    ConstructorArgs&&...)>;

    template <typename ShapeFunction>
    using IntegrationMethod = typename NumLib::GaussLegendreIntegrationPolicy<
        typename ShapeFunction::MeshElement>::IntegrationMethod;

    template <typename ShapeFunctionDisplacement,
              typename ShapeFunctionPressure>
    using LAData =
        LocalAssemblerData<ShapeFunctionDisplacement, ShapeFunctionPressure,
                           IntegrationMethod<ShapeFunctionDisplacement>,
                           GlobalDim>;

    /// Generates a function that creates a new LocalAssembler of type
    /// LAData<ShapeFunctionDisplacement, ShapeFunctionPressure>. Only functions
    /// with shape function's dimension less or equal to the global dimension
    /// are instantiated, e.g.  following combinations of shape functions and
    /// global dimensions: (Line2, 1), (Line2, 2), (Line2, 3), (Hex20, 3) but
    /// not (Hex20, 2) or (Hex20, 1).
    template <typename ShapeFunction, typename LowerOrderShapeFunction>
    static LADataBuilder makeLocalAssemblerBuilder()
    {
        return [](MeshLib::Element const& e,
                  std::size_t const local_matrix_size,
                  ConstructorArgs&&... args)
        {
            return LADataIntfPtr{
                new LAData<ShapeFunction, LowerOrderShapeFunction>{
                    e, local_matrix_size,
                    std::forward<ConstructorArgs>(args)...}};
        };
    }

    /// Mapping of element types to local assembler constructors.
    std::unordered_map<std::type_index, LADataBuilder> _builder;

    NumLib::LocalToGlobalIndexMap const& _dof_table;
};

}  // namespace ProcessLib::StokesFlow
