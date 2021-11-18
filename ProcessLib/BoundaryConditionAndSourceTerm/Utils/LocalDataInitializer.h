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
namespace BoundaryConditionAndSourceTerm
{
/// The LocalDataInitializer is a functor creating a local assembler data with
/// corresponding to the mesh element type shape functions and calling
/// initialization of the new local assembler data.
/// For example for MeshLib::Quad a local assembler data with template argument
/// NumLib::ShapeQuad4 is created.
template <typename LocalAssemblerInterface,
          template <typename, typename, int> class LocalAssemblerData,
          int GlobalDim, typename... ConstructorArgs>
class LocalDataInitializer final
{
    struct Is2ndOrderElement
    {
        template <typename ElementTraits>
        constexpr bool operator()(ElementTraits*) const
        {
            return ElementTraits::ShapeFunction::ORDER == 2 ||
                   // points are needed for 2nd order, too
                   std::is_same_v<MeshLib::Point,
                                  typename ElementTraits::Element>;
        }
    };

public:
    using LADataIntfPtr = std::unique_ptr<LocalAssemblerInterface>;

    LocalDataInitializer(NumLib::LocalToGlobalIndexMap const& dof_table,
                         const unsigned shapefunction_order)
        : _dof_table(dof_table)
    {
        if (shapefunction_order < 1 || 2 < shapefunction_order)
            OGS_FATAL("The given shape function order {:d} is not supported",
                      shapefunction_order);

        if (shapefunction_order == 1)
        {
            // 1st order is enabled on all elements
            using EnabledElementTraits = EnabledElementTraitsLagrange;

            BaseLib::TMP::foreach<EnabledElementTraits>(
                [this]<typename ET>(ET*)
                {
                    using Elt = typename ET::Element;
                    // this will use linear shape functions on higher order
                    // elements and the linear shape function on linear elements
                    using LowShp = typename ET::LowerOrderShapeFunction;
                    _builder[std::type_index(typeid(Elt))] =
                        makeLocalAssemblerBuilder<LowShp>();
                });
        }
        else if (shapefunction_order == 2)
        {
            // 2nd order only on 2nd order elements
            using EnabledElementTraits =
                decltype(BaseLib::TMP::filter<EnabledElementTraitsLagrange>(
                    std::declval<Is2ndOrderElement>()));

            BaseLib::TMP::foreach<EnabledElementTraits>(
                [this]<typename ET>(ET*)
                {
                    using Elt = typename ET::Element;
                    using Shp2ndOrder = typename ET::ShapeFunction;
                    _builder[std::type_index(typeid(Elt))] =
                        makeLocalAssemblerBuilder<Shp2ndOrder>();
                });
        }
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

    template <typename ShapeFunction>
    using LAData =
        LocalAssemblerData<ShapeFunction, IntegrationMethod<ShapeFunction>,
                           GlobalDim>;

    /// A helper forwarding to the correct version of makeLocalAssemblerBuilder
    /// depending whether the global dimension is less than the shape function's
    /// dimension or not.
    template <typename ShapeFunction>
    static LADataBuilder makeLocalAssemblerBuilder()
    {
        return makeLocalAssemblerBuilder<ShapeFunction>(
            static_cast<std::integral_constant<
                bool, (GlobalDim >= ShapeFunction::DIM)>*>(nullptr));
    }

    /// Mapping of element types to local assembler constructors.
    std::unordered_map<std::type_index, LADataBuilder> _builder;

    NumLib::LocalToGlobalIndexMap const& _dof_table;

    // local assembler builder implementations.
private:
    /// Generates a function that creates a new LocalAssembler of type
    /// LAData<ShapeFunction>. Only functions with shape function's dimension
    /// less or equal to the global dimension are instantiated, e.g.  following
    /// combinations of shape functions and global dimensions: (Line2, 1),
    /// (Line2, 2), (Line2, 3), (Hex20, 3) but not (Hex20, 2) or (Hex20, 1).
    template <typename ShapeFunction>
    static LADataBuilder makeLocalAssemblerBuilder(std::true_type* /*unused*/)
    {
        return [](MeshLib::Element const& e,
                  std::size_t const local_matrix_size,
                  ConstructorArgs&&... args)
        {
            return LADataIntfPtr{new LAData<ShapeFunction>{
                e, local_matrix_size, std::forward<ConstructorArgs>(args)...}};
        };
    }

    /// Returns nullptr for shape functions whose dimensions are less than the
    /// global dimension.
    template <typename ShapeFunction>
    static LADataBuilder makeLocalAssemblerBuilder(std::false_type* /*unused*/)
    {
        return nullptr;
    }
};

}  // namespace BoundaryConditionAndSourceTerm
}  // namespace ProcessLib
