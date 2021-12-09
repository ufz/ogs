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
#include "ProcessLib/HeatTransportBHE/BHE/BHETypes.h"
#include "ProcessLib/Utils/EnabledElements.h"

namespace ProcessLib
{
namespace HeatTransportBHE
{
/// The LocalDataInitializer is a functor creating a local assembler data with
/// corresponding to the mesh element type shape functions and calling
/// initialization of the new local assembler data.
/// For example for MeshLib::Quad a local assembler data with template argument
/// NumLib::ShapeQuad4 is created.
template <typename LocalAssemblerInterface,
          template <typename, typename>
          class LocalAssemblerDataSoil,
          template <typename, typename, typename>
          class LocalAssemblerDataBHE,
          typename... ConstructorArgs>
class LocalDataInitializer final
{
    template <unsigned DIM>
    struct IsNDElement
    {
        template <typename ElementTraits>
        constexpr bool operator()(ElementTraits*) const
        {
            return ElementTraits::Element::dimension == DIM;
        }
    };

public:
    using LADataIntfPtr = std::unique_ptr<LocalAssemblerInterface>;

    explicit LocalDataInitializer(
        NumLib::LocalToGlobalIndexMap const& dof_table)
        : _dof_table(dof_table)
    {
        // 3D soil elements
        using Enabled3DElementTraits =
            decltype(BaseLib::TMP::filter<EnabledElementTraitsLagrange>(
                std::declval<IsNDElement<3>>()));

        BaseLib::TMP::foreach<Enabled3DElementTraits>(
            [this]<typename ET>(ET*)
            {
                using MeshElement = typename ET::Element;
                using ShapeFunction = typename ET::ShapeFunction;

                _builder[std::type_index(typeid(MeshElement))] =
                    makeLocalAssemblerBuilder<ShapeFunction>();
            });

        // 1D BHE elements
        using Enabled1DElementTraits =
            decltype(BaseLib::TMP::filter<EnabledElementTraitsLagrange>(
                std::declval<IsNDElement<1>>()));

        BaseLib::TMP::foreach<Enabled1DElementTraits>(
            [this]<typename ET>(ET*)
            {
                using MeshElement = typename ET::Element;
                using ShapeFunction = typename ET::ShapeFunction;

                _builder[std::type_index(typeid(MeshElement))] =
                    makeLocalAssemblerBuilderBHE<ShapeFunction>();
            });
    }

    /// Returns data pointer to the newly created local assembler data.
    ///
    /// \attention
    /// The index \c id is not necessarily the mesh item's id. Especially when
    /// having multiple meshes it will differ from the latter.
    LADataIntfPtr operator()(
        std::size_t const /*id*/,
        MeshLib::Element const& mesh_item,
        std::unordered_map<std::size_t, BHE::BHETypes*> const&
            element_to_bhe_map,
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

        return it->second(mesh_item,
                          element_to_bhe_map,
                          std::forward<ConstructorArgs>(args)...);
    }

private:
    using LADataBuilder = std::function<LADataIntfPtr(
        MeshLib::Element const& e,
        std::unordered_map<std::size_t, BHE::BHETypes*> const&
            element_to_bhe_map,
        ConstructorArgs&&...)>;

    template <typename ShapeFunction>
    using IntegrationMethod = typename NumLib::GaussLegendreIntegrationPolicy<
        typename ShapeFunction::MeshElement>::IntegrationMethod;

    // local assembler builder implementations.
    template <typename ShapeFunction>
    using LADataSoil =
        LocalAssemblerDataSoil<ShapeFunction, IntegrationMethod<ShapeFunction>>;

    template <typename ShapeFunction>
    static LADataBuilder makeLocalAssemblerBuilder()
    {
        return [](MeshLib::Element const& e,
                  std::unordered_map<std::size_t, BHE::BHETypes*> const&
                  /* unused */,
                  ConstructorArgs&&... args) -> LADataIntfPtr {
            if (e.getDimension() == 3)  // soil elements
            {
                return LADataIntfPtr{new LADataSoil<ShapeFunction>{
                    e, std::forward<ConstructorArgs>(args)...}};
            }

            return nullptr;
        };
    }

    template <typename ShapeFunction, typename BHEType>
    using LADataBHE = LocalAssemblerDataBHE<ShapeFunction,
                                            IntegrationMethod<ShapeFunction>,
                                            BHEType>;
    template <typename ShapeFunction>
    static LADataBuilder makeLocalAssemblerBuilderBHE()
    {
        return [](MeshLib::Element const& e,
                  std::unordered_map<std::size_t, BHE::BHETypes*> const&
                      element_to_bhe_map,
                  ConstructorArgs&&... args) -> LADataIntfPtr {
            auto& bhe = *element_to_bhe_map.at(e.getID());

            if (std::holds_alternative<BHE::BHE_1U>(bhe))
            {
                return LADataIntfPtr{new LADataBHE<ShapeFunction, BHE::BHE_1U>{
                    e, std::get<BHE::BHE_1U>(bhe),
                    std::forward<ConstructorArgs>(args)...}};
            }

            if (std::holds_alternative<BHE::BHE_CXA>(bhe))
            {
                return LADataIntfPtr{new LADataBHE<ShapeFunction, BHE::BHE_CXA>{
                    e, std::get<BHE::BHE_CXA>(bhe),
                    std::forward<ConstructorArgs>(args)...}};
            }

            if (std::holds_alternative<BHE::BHE_CXC>(bhe))
            {
                return LADataIntfPtr{new LADataBHE<ShapeFunction, BHE::BHE_CXC>{
                    e, std::get<BHE::BHE_CXC>(bhe),
                    std::forward<ConstructorArgs>(args)...}};
            }

            if (std::holds_alternative<BHE::BHE_2U>(bhe))
            {
                return LADataIntfPtr{new LADataBHE<ShapeFunction, BHE::BHE_2U>{
                    e, std::get<BHE::BHE_2U>(bhe),
                    std::forward<ConstructorArgs>(args)...}};
            }

            if (std::holds_alternative<BHE::BHE_1P>(bhe))
            {
                return LADataIntfPtr{new LADataBHE<ShapeFunction, BHE::BHE_1P>{
                    e, std::get<BHE::BHE_1P>(bhe),
                    std::forward<ConstructorArgs>(args)...}};
            }
            OGS_FATAL(
                "Trying to create local assembler for an unknown BHE type.");
        };
    }

    /// Mapping of element types to local assembler constructors.
    std::unordered_map<std::type_index, LADataBuilder> _builder;

    NumLib::LocalToGlobalIndexMap const& _dof_table;
};
}  // namespace HeatTransportBHE
}  // namespace ProcessLib
