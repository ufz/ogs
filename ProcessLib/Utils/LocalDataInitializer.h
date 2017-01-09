/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESS_LIB_LOCALDATAINITIALIVER_H_
#define PROCESS_LIB_LOCALDATAINITIALIVER_H_

#include <functional>
#include <memory>
#include <typeindex>
#include <typeinfo>
#include <type_traits>
#include <unordered_map>

#include "MeshLib/Elements/Elements.h"
#include "NumLib/DOF/LocalToGlobalIndexMap.h"
#include "NumLib/Fem/Integration/GaussIntegrationPolicy.h"

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
// Faces of tets, pyramids and prisms are triangles
#define ENABLED_ELEMENT_TYPE_TRI                                       \
    ((ENABLED_ELEMENT_TYPE_SIMPLEX) | (ENABLED_ELEMENT_TYPE_PYRAMID) | \
     (ENABLED_ELEMENT_TYPE_PRISM))
// Faces of hexes, pyramids and prisms are quads
#define ENABLED_ELEMENT_TYPE_QUAD                                     \
    ((ENABLED_ELEMENT_TYPE_CUBOID) | (ENABLED_ELEMENT_TYPE_PYRAMID) | \
     (ENABLED_ELEMENT_TYPE_PRISM))
// Faces of triangles and quads are lines
#define ENABLED_ELEMENT_TYPE_LINE \
    ((ENABLED_ELEMENT_TYPE_TRI) | (ENABLED_ELEMENT_TYPE_QUAD))

// All enabled element types
#define OGS_ENABLED_ELEMENTS                                          \
    ((ENABLED_ELEMENT_TYPE_SIMPLEX) | (ENABLED_ELEMENT_TYPE_CUBOID) | \
     (ENABLED_ELEMENT_TYPE_PYRAMID) | (ENABLED_ELEMENT_TYPE_PRISM))

// Include only what is needed (Well, the conditions are not sharp).
#if (OGS_ENABLED_ELEMENTS & ENABLED_ELEMENT_TYPE_LINE) != 0
#include "NumLib/Fem/ShapeFunction/ShapeLine2.h"
#include "NumLib/Fem/ShapeFunction/ShapeLine3.h"
#include "NumLib/Fem/ShapeFunction/ShapePoint1.h"
#endif

#if (OGS_ENABLED_ELEMENTS & ENABLED_ELEMENT_TYPE_SIMPLEX) != 0
#include "NumLib/Fem/ShapeFunction/ShapeTet10.h"
#include "NumLib/Fem/ShapeFunction/ShapeTet4.h"
#endif

#if (OGS_ENABLED_ELEMENTS & ENABLED_ELEMENT_TYPE_TRI) != 0
#include "NumLib/Fem/ShapeFunction/ShapeTri3.h"
#include "NumLib/Fem/ShapeFunction/ShapeTri6.h"
#endif

#if (OGS_ENABLED_ELEMENTS & ENABLED_ELEMENT_TYPE_CUBOID) != 0
#include "NumLib/Fem/ShapeFunction/ShapeHex20.h"
#include "NumLib/Fem/ShapeFunction/ShapeHex8.h"
#endif

#if (OGS_ENABLED_ELEMENTS & ENABLED_ELEMENT_TYPE_QUAD) != 0
#include "NumLib/Fem/ShapeFunction/ShapeQuad4.h"
#include "NumLib/Fem/ShapeFunction/ShapeQuad8.h"
#include "NumLib/Fem/ShapeFunction/ShapeQuad9.h"
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
/// The LocalDataInitializer is a functor creating a local assembler data with
/// corresponding to the mesh element type shape functions and calling
/// initialization of the new local assembler data.
/// For example for MeshLib::Line a local assembler data with template argument
/// NumLib::ShapeLine2 is created.
template <typename LocalAssemblerInterface,
          template <typename, typename, unsigned> class LocalAssemblerData,
          unsigned GlobalDim, typename... ConstructorArgs>
class LocalDataInitializer final
{
public:
    using LADataIntfPtr = std::unique_ptr<LocalAssemblerInterface>;

    LocalDataInitializer(NumLib::LocalToGlobalIndexMap const& dof_table,
                         const unsigned shapefunction_order)
        : _dof_table(dof_table)
    {
        if (shapefunction_order < 1 || 2 < shapefunction_order)
            OGS_FATAL("The given shape function order %d is not supported",
                      shapefunction_order);

        if (shapefunction_order == 1)
        {
// /// Lines and points ///////////////////////////////////

#if (OGS_ENABLED_ELEMENTS & ENABLED_ELEMENT_TYPE_LINE) != 0 && \
    OGS_MAX_ELEMENT_DIM >= 0 && OGS_MAX_ELEMENT_ORDER >= 1
            _builder[std::type_index(typeid(MeshLib::Point))] =
                makeLocalAssemblerBuilder<NumLib::ShapePoint1>();
#endif

#if (OGS_ENABLED_ELEMENTS & ENABLED_ELEMENT_TYPE_LINE) != 0 && \
    OGS_MAX_ELEMENT_DIM >= 1 && OGS_MAX_ELEMENT_ORDER >= 1
            _builder[std::type_index(typeid(MeshLib::Line))] =
                makeLocalAssemblerBuilder<NumLib::ShapeLine2>();
#endif

#if (OGS_ENABLED_ELEMENTS & ENABLED_ELEMENT_TYPE_LINE) != 0 && \
    OGS_MAX_ELEMENT_DIM >= 1 && OGS_MAX_ELEMENT_ORDER >= 2
            _builder[std::type_index(typeid(MeshLib::Line3))] =
                makeLocalAssemblerBuilder<NumLib::ShapeLine2>();
#endif

// /// Quads and Hexahedra ///////////////////////////////////

#if (OGS_ENABLED_ELEMENTS & ENABLED_ELEMENT_TYPE_QUAD) != 0 && \
    OGS_MAX_ELEMENT_DIM >= 2 && OGS_MAX_ELEMENT_ORDER >= 1
            _builder[std::type_index(typeid(MeshLib::Quad))] =
                makeLocalAssemblerBuilder<NumLib::ShapeQuad4>();
#endif

#if (OGS_ENABLED_ELEMENTS & ENABLED_ELEMENT_TYPE_CUBOID) != 0 && \
    OGS_MAX_ELEMENT_DIM >= 3 && OGS_MAX_ELEMENT_ORDER >= 1
            _builder[std::type_index(typeid(MeshLib::Hex))] =
                makeLocalAssemblerBuilder<NumLib::ShapeHex8>();
#endif

#if (OGS_ENABLED_ELEMENTS & ENABLED_ELEMENT_TYPE_QUAD) != 0 && \
    OGS_MAX_ELEMENT_DIM >= 2 && OGS_MAX_ELEMENT_ORDER >= 2
            _builder[std::type_index(typeid(MeshLib::Quad8))] =
                makeLocalAssemblerBuilder<NumLib::ShapeQuad4>();
            _builder[std::type_index(typeid(MeshLib::Quad9))] =
                makeLocalAssemblerBuilder<NumLib::ShapeQuad4>();
#endif

#if (OGS_ENABLED_ELEMENTS & ENABLED_ELEMENT_TYPE_CUBOID) != 0 && \
    OGS_MAX_ELEMENT_DIM >= 3 && OGS_MAX_ELEMENT_ORDER >= 2
            _builder[std::type_index(typeid(MeshLib::Hex20))] =
                makeLocalAssemblerBuilder<NumLib::ShapeHex8>();
#endif

// /// Simplices ////////////////////////////////////////////////

#if (OGS_ENABLED_ELEMENTS & ENABLED_ELEMENT_TYPE_TRI) != 0 && \
    OGS_MAX_ELEMENT_DIM >= 2 && OGS_MAX_ELEMENT_ORDER >= 1
            _builder[std::type_index(typeid(MeshLib::Tri))] =
                makeLocalAssemblerBuilder<NumLib::ShapeTri3>();
#endif

#if (OGS_ENABLED_ELEMENTS & ENABLED_ELEMENT_TYPE_SIMPLEX) != 0 && \
    OGS_MAX_ELEMENT_DIM >= 3 && OGS_MAX_ELEMENT_ORDER >= 1
            _builder[std::type_index(typeid(MeshLib::Tet))] =
                makeLocalAssemblerBuilder<NumLib::ShapeTet4>();
#endif

#if (OGS_ENABLED_ELEMENTS & ENABLED_ELEMENT_TYPE_TRI) != 0 && \
    OGS_MAX_ELEMENT_DIM >= 2 && OGS_MAX_ELEMENT_ORDER >= 2
            _builder[std::type_index(typeid(MeshLib::Tri6))] =
                makeLocalAssemblerBuilder<NumLib::ShapeTri3>();
#endif

#if (OGS_ENABLED_ELEMENTS & ENABLED_ELEMENT_TYPE_SIMPLEX) != 0 && \
    OGS_MAX_ELEMENT_DIM >= 3 && OGS_MAX_ELEMENT_ORDER >= 2
            _builder[std::type_index(typeid(MeshLib::Tet10))] =
                makeLocalAssemblerBuilder<NumLib::ShapeTet4>();
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
                makeLocalAssemblerBuilder<NumLib::ShapePrism6>();
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
                makeLocalAssemblerBuilder<NumLib::ShapePyra5>();
#endif
        }
        else if (shapefunction_order == 2)
        {
// /// Lines and points ///////////////////////////////////

#if (OGS_ENABLED_ELEMENTS & ENABLED_ELEMENT_TYPE_LINE) != 0 && \
    OGS_MAX_ELEMENT_DIM >= 1 && OGS_MAX_ELEMENT_ORDER >= 2
            _builder[std::type_index(typeid(MeshLib::Line3))] =
                makeLocalAssemblerBuilder<NumLib::ShapeLine3>();
#endif

// /// Quads and Hexahedra ///////////////////////////////////

#if (OGS_ENABLED_ELEMENTS & ENABLED_ELEMENT_TYPE_QUAD) != 0 && \
    OGS_MAX_ELEMENT_DIM >= 2 && OGS_MAX_ELEMENT_ORDER >= 2
            _builder[std::type_index(typeid(MeshLib::Quad8))] =
                makeLocalAssemblerBuilder<NumLib::ShapeQuad8>();
            _builder[std::type_index(typeid(MeshLib::Quad9))] =
                makeLocalAssemblerBuilder<NumLib::ShapeQuad9>();
#endif

#if (OGS_ENABLED_ELEMENTS & ENABLED_ELEMENT_TYPE_CUBOID) != 0 && \
    OGS_MAX_ELEMENT_DIM >= 3 && OGS_MAX_ELEMENT_ORDER >= 2
            _builder[std::type_index(typeid(MeshLib::Hex20))] =
                makeLocalAssemblerBuilder<NumLib::ShapeHex20>();
#endif

// /// Simplices ////////////////////////////////////////////////

#if (OGS_ENABLED_ELEMENTS & ENABLED_ELEMENT_TYPE_TRI) != 0 && \
    OGS_MAX_ELEMENT_DIM >= 2 && OGS_MAX_ELEMENT_ORDER >= 2
            _builder[std::type_index(typeid(MeshLib::Tri6))] =
                makeLocalAssemblerBuilder<NumLib::ShapeTri6>();
#endif

#if (OGS_ENABLED_ELEMENTS & ENABLED_ELEMENT_TYPE_SIMPLEX) != 0 && \
    OGS_MAX_ELEMENT_DIM >= 3 && OGS_MAX_ELEMENT_ORDER >= 2
            _builder[std::type_index(typeid(MeshLib::Tet10))] =
                makeLocalAssemblerBuilder<NumLib::ShapeTet10>();
#endif

// /// Prisms ////////////////////////////////////////////////////

#if (OGS_ENABLED_ELEMENTS & ENABLED_ELEMENT_TYPE_PRISM) != 0 && \
    OGS_MAX_ELEMENT_DIM >= 3 && OGS_MAX_ELEMENT_ORDER >= 2
            _builder[std::type_index(typeid(MeshLib::Prism15))] =
                makeLocalAssemblerBuilder<NumLib::ShapePrism15>();
#endif

// /// Pyramids //////////////////////////////////////////////////

#if (OGS_ENABLED_ELEMENTS & ENABLED_ELEMENT_TYPE_PYRAMID) != 0 && \
    OGS_MAX_ELEMENT_DIM >= 3 && OGS_MAX_ELEMENT_ORDER >= 2
            _builder[std::type_index(typeid(MeshLib::Pyramid13))] =
                makeLocalAssemblerBuilder<NumLib::ShapePyra13>();
#endif
        }
    }

    /// Sets the provided \c data_ptr to the newly created local assembler data.
    ///
    /// \attention
    /// The index \c id is not necessarily the mesh item's id. Especially when
    /// having multiple meshes it will differ from the latter.
    void operator()(std::size_t const id,
                    MeshLib::Element const& mesh_item,
                    LADataIntfPtr& data_ptr,
                    ConstructorArgs&&... args) const
    {
        auto const type_idx = std::type_index(typeid(mesh_item));
        auto const it = _builder.find(type_idx);

        if (it != _builder.end())
        {
            auto const num_local_dof = _dof_table.getNumberOfElementDOF(id);
            data_ptr = it->second(mesh_item, num_local_dof,
                                  std::forward<ConstructorArgs>(args)...);
        }
        else
        {
            OGS_FATAL(
                "You are trying to build a local assembler for an unknown mesh "
                "element type (%s)."
                " Maybe you have disabled this mesh element type in your build "
                "configuration.",
                type_idx.name());
        }
    }

private:
    using LADataBuilder =
        std::function<LADataIntfPtr(MeshLib::Element const& e,
                                    std::size_t const local_matrix_size,
                                    ConstructorArgs&&...)>;

    template <typename ShapeFunction>
    using IntegrationMethod = typename NumLib::GaussIntegrationPolicy<
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
    static LADataBuilder makeLocalAssemblerBuilder(std::true_type*)
    {
        return [](MeshLib::Element const& e,
                  std::size_t const local_matrix_size,
                  ConstructorArgs&&... args) {
            return LADataIntfPtr{new LAData<ShapeFunction>{
                e, local_matrix_size, std::forward<ConstructorArgs>(args)...}};
        };
    }

    /// Returns nullptr for shape functions whose dimensions are less than the
    /// global dimension.
    template <typename ShapeFunction>
    static LADataBuilder makeLocalAssemblerBuilder(std::false_type*)
    {
        return nullptr;
    }
};

}  // namespace ProcessLib

#undef ENABLED_ELEMENT_TYPE_SIMPLEX
#undef ENABLED_ELEMENT_TYPE_CUBOID
#undef ENABLED_ELEMENT_TYPE_PYRAMID
#undef ENABLED_ELEMENT_TYPE_PRISM
#undef ENABLED_ELEMENT_TYPE_LINE
#undef ENABLED_ELEMENT_TYPE_TRI
#undef ENABLED_ELEMENT_TYPE_QUAD
#undef OGS_ENABLED_ELEMENTS

#endif  // PROCESS_LIB_LOCALDATAINITIALIVER_H_
