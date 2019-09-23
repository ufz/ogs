/**
 * \file
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
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

#ifndef OGS_MAX_ELEMENT_DIM
static_assert(false, "The macro OGS_MAX_ELEMENT_DIM is undefined.");
#endif

#ifndef OGS_MAX_ELEMENT_ORDER
static_assert(false, "The macro OGS_MAX_ELEMENT_ORDER is undefined.");
#endif

// The following macros decide which element types will be compiled, i.e.
// which element types will be available for use in simulations.

#ifdef OGS_ENABLE_ELEMENT_SIMPLEX
#define ENABLED_ELEMENT_TYPE_SIMPLEX 1u
#else
#define ENABLED_ELEMENT_TYPE_SIMPLEX 0u
#endif

#ifdef OGS_ENABLE_ELEMENT_CUBOID
#define ENABLED_ELEMENT_TYPE_CUBOID 1u << 1
#else
#define ENABLED_ELEMENT_TYPE_CUBOID 0u
#endif

#ifdef OGS_ENABLE_ELEMENT_PRISM
#define ENABLED_ELEMENT_TYPE_PRISM 1u << 2
#else
#define ENABLED_ELEMENT_TYPE_PRISM 0u
#endif

#ifdef OGS_ENABLE_ELEMENT_PYRAMID
#define ENABLED_ELEMENT_TYPE_PYRAMID 1u << 3
#else
#define ENABLED_ELEMENT_TYPE_PYRAMID 0u
#endif

// Dependent element types.
// All enabled element types
#define OGS_ENABLED_ELEMENTS                                          \
    ((ENABLED_ELEMENT_TYPE_SIMPLEX) | (ENABLED_ELEMENT_TYPE_CUBOID) | \
     (ENABLED_ELEMENT_TYPE_PYRAMID) | (ENABLED_ELEMENT_TYPE_PRISM))

// Include only what is needed (Well, the conditions are not sharp).
#if OGS_ENABLED_ELEMENTS != 0
#include "NumLib/Fem/ShapeFunction/ShapeLine2.h"
#include "NumLib/Fem/ShapeFunction/ShapeLine3.h"
#endif

#if (OGS_ENABLED_ELEMENTS & ENABLED_ELEMENT_TYPE_SIMPLEX) != 0
#include "NumLib/Fem/ShapeFunction/ShapeTet10.h"
#include "NumLib/Fem/ShapeFunction/ShapeTet4.h"
#endif

#if (OGS_ENABLED_ELEMENTS & ENABLED_ELEMENT_TYPE_CUBOID) != 0
#include "NumLib/Fem/ShapeFunction/ShapeHex20.h"
#include "NumLib/Fem/ShapeFunction/ShapeHex8.h"
#endif

#if (OGS_ENABLED_ELEMENTS & ENABLED_ELEMENT_TYPE_PRISM) != 0
#include "NumLib/Fem/ShapeFunction/ShapePrism15.h"
#include "NumLib/Fem/ShapeFunction/ShapePrism6.h"
#endif

#if (OGS_ENABLED_ELEMENTS & ENABLED_ELEMENT_TYPE_PYRAMID) != 0
#include "NumLib/Fem/ShapeFunction/ShapePyra13.h"
#include "NumLib/Fem/ShapeFunction/ShapePyra5.h"
#endif

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
public:
    using LADataIntfPtr = std::unique_ptr<LocalAssemblerInterface>;

    explicit LocalDataInitializer(
        NumLib::LocalToGlobalIndexMap const& dof_table)
        : _dof_table(dof_table)
    {
        // REMARKS: At the moment, only a 3D mesh (soil) with 1D elements (BHE)
        // are supported.
#if (OGS_ENABLED_ELEMENTS & ENABLED_ELEMENT_TYPE_CUBOID) != 0 && \
    OGS_MAX_ELEMENT_DIM >= 3 && OGS_MAX_ELEMENT_ORDER >= 1
        _builder[std::type_index(typeid(MeshLib::Hex))] =
            makeLocalAssemblerBuilder<NumLib::ShapeHex8>();
#endif

#if (OGS_ENABLED_ELEMENTS & ENABLED_ELEMENT_TYPE_CUBOID) != 0 && \
    OGS_MAX_ELEMENT_DIM >= 3 && OGS_MAX_ELEMENT_ORDER >= 2
        _builder[std::type_index(typeid(MeshLib::Hex20))] =
            makeLocalAssemblerBuilder<NumLib::ShapeHex20>();
#endif

        // /// Simplices ////////////////////////////////////////////////
#if (OGS_ENABLED_ELEMENTS & ENABLED_ELEMENT_TYPE_SIMPLEX) != 0 && \
    OGS_MAX_ELEMENT_DIM >= 3 && OGS_MAX_ELEMENT_ORDER >= 1
        _builder[std::type_index(typeid(MeshLib::Tet))] =
            makeLocalAssemblerBuilder<NumLib::ShapeTet4>();
#endif

#if (OGS_ENABLED_ELEMENTS & ENABLED_ELEMENT_TYPE_SIMPLEX) != 0 && \
    OGS_MAX_ELEMENT_DIM >= 3 && OGS_MAX_ELEMENT_ORDER >= 2
        _builder[std::type_index(typeid(MeshLib::Tet10))] =
            makeLocalAssemblerBuilder<NumLib::ShapeTet10>();
#endif

        // /// Prisms ////////////////////////////////////////////////////

#if (OGS_ENABLED_ELEMENTS & ENABLED_ELEMENT_TYPE_PRISM) != 0 && \
    OGS_MAX_ELEMENT_DIM >= 3 && OGS_MAX_ELEMENT_ORDER >= 1
        _builder[std::type_index(typeid(MeshLib::Prism))] =
            makeLocalAssemblerBuilder<NumLib::ShapePrism6>();
#endif

#if (OGS_ENABLED_ELEMENTS & ENABLED_ELEMENT_TYPE_PRISM) != 0 && \
    OGS_MAX_ELEMENT_DIM >= 3 && OGS_MAX_ELEMENT_ORDER >= 2
        _builder[std::type_index(typeid(MeshLib::Prism15))] =
            makeLocalAssemblerBuilder<NumLib::ShapePrism15>();
#endif

        // /// Pyramids //////////////////////////////////////////////////

#if (OGS_ENABLED_ELEMENTS & ENABLED_ELEMENT_TYPE_PYRAMID) != 0 && \
    OGS_MAX_ELEMENT_DIM >= 3 && OGS_MAX_ELEMENT_ORDER >= 1
        _builder[std::type_index(typeid(MeshLib::Pyramid))] =
            makeLocalAssemblerBuilder<NumLib::ShapePyra5>();
#endif

#if (OGS_ENABLED_ELEMENTS & ENABLED_ELEMENT_TYPE_PYRAMID) != 0 && \
    OGS_MAX_ELEMENT_DIM >= 3 && OGS_MAX_ELEMENT_ORDER >= 2
        _builder[std::type_index(typeid(MeshLib::Pyramid13))] =
            makeLocalAssemblerBuilder<NumLib::ShapePyra13>();
#endif
        // /// Lines ///////////////////////////////////

#if OGS_MAX_ELEMENT_DIM >= 2 && OGS_MAX_ELEMENT_ORDER >= 1
        _builder[std::type_index(typeid(MeshLib::Line))] =
            makeLocalAssemblerBuilderBHE<NumLib::ShapeLine2>();
#endif

#if OGS_MAX_ELEMENT_DIM >= 3 && OGS_MAX_ELEMENT_ORDER >= 2
        _builder[std::type_index(typeid(MeshLib::Line3))] =
            makeLocalAssemblerBuilderBHE<NumLib::ShapeLine3>();
#endif
    }

    /// Sets the provided \c data_ptr to the newly created local assembler data.
    ///
    /// \attention
    /// The index \c id is not necessarily the mesh item's id. Especially when
    /// having multiple meshes it will differ from the latter.
    void operator()(std::size_t const /*id*/,
                    MeshLib::Element const& mesh_item,
                    LADataIntfPtr& data_ptr,
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
                "element type (%s)."
                " Maybe you have disabled this mesh element type in your build "
                "configuration or this process requires higher order elements.",
                type_idx.name());
        }

        data_ptr = it->second(mesh_item,
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
            OGS_FATAL(
                "Trying to create local assembler for an unknown BHE type.");
        };
    }

    /// Mapping of element types to local assembler constructors.
    std::unordered_map<std::type_index, LADataBuilder> _builder;

    NumLib::LocalToGlobalIndexMap const& _dof_table;
};  // namespace HeatTransportBHE
}  // namespace HeatTransportBHE
}  // namespace ProcessLib

#undef ENABLED_ELEMENT_TYPE_SIMPLEX
#undef ENABLED_ELEMENT_TYPE_CUBOID
#undef ENABLED_ELEMENT_TYPE_PYRAMID
#undef ENABLED_ELEMENT_TYPE_PRISM
#undef OGS_ENABLED_ELEMENTS
