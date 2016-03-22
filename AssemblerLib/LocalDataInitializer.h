/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef ASSEMBLER_LIB_LOCALDATAINITIALIZER_H_
#define ASSEMBLER_LIB_LOCALDATAINITIALIZER_H_

#include <functional>
#include <typeindex>
#include <typeinfo>
#include <unordered_map>

#include "MeshLib/Elements/Elements.h"

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
#define ENABLED_ELEMENT_TYPE_PRISM   0u
#endif

#ifdef OGS_ENABLE_ELEMENT_PYRAMID
#define ENABLED_ELEMENT_TYPE_PYRAMID 1u << 3
#else
#define ENABLED_ELEMENT_TYPE_PYRAMID 0u
#endif

// Dependent element types.
// Faces of tets, pyramids and prisms are triangles
#define ENABLED_ELEMENT_TYPE_TRI  ((ENABLED_ELEMENT_TYPE_SIMPLEX) | (ENABLED_ELEMENT_TYPE_PYRAMID) \
                                  |(ENABLED_ELEMENT_TYPE_PRISM))
// Faces of hexes, pyramids and prisms are quads
#define ENABLED_ELEMENT_TYPE_QUAD ((ENABLED_ELEMENT_TYPE_CUBOID) | (ENABLED_ELEMENT_TYPE_PYRAMID) \
                                  |(ENABLED_ELEMENT_TYPE_PRISM))
// Faces of triangles and quads are lines
#define ENABLED_ELEMENT_TYPE_LINE ((ENABLED_ELEMENT_TYPE_TRI) | (ENABLED_ELEMENT_TYPE_QUAD))

// All enabled element types
#define OGS_ENABLED_ELEMENTS \
    ( (ENABLED_ELEMENT_TYPE_SIMPLEX) \
    | (ENABLED_ELEMENT_TYPE_CUBOID ) \
    | (ENABLED_ELEMENT_TYPE_PYRAMID) \
    | (ENABLED_ELEMENT_TYPE_PRISM  ) )


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


// Generates a lambda that creates a new LocalAssembler of type LAData<SHAPE_FCT>
#define OGS_NEW_LOCAL_ASSEMBLER(SHAPE_FCT) \
    [](MeshLib::Element const& e, \
       std::size_t const local_matrix_size, \
       unsigned const integration_order, \
       ConstructorArgs&&... args) \
    { \
        return new LAData<SHAPE_FCT>{ \
            e, local_matrix_size, integration_order, \
            std::forward<ConstructorArgs>(args)... \
        }; \
    }


namespace AssemblerLib
{

/// The LocalDataInitializer is a functor creating a local assembler data with
/// corresponding to the mesh element type shape functions and calling
/// initialization of the new local assembler data.
/// For example for MeshLib::Line a local assembler data with template argument
/// NumLib::ShapeLine2 is created.
template <
    template <typename, typename> class LocalAssemblerDataInterface_,
    template <typename, typename, typename, typename, unsigned> class LocalAssemblerData_,
    typename GlobalMatrix_,
    typename GlobalVector_,
    unsigned GlobalDim,
    typename... ConstructorArgs>
class LocalDataInitializer
{
    template <typename ShapeFunction_>
    using IntegrationMethod = typename NumLib::GaussIntegrationPolicy<
                typename ShapeFunction_::MeshElement>::IntegrationMethod;

    template <typename ShapeFunction_>
        using LAData = LocalAssemblerData_<
                ShapeFunction_,
                IntegrationMethod<ShapeFunction_>,
                GlobalMatrix_, GlobalVector_, GlobalDim>;

public:
    LocalDataInitializer()
    {
        // /// Lines and points ///////////////////////////////////

#if (OGS_ENABLED_ELEMENTS & ENABLED_ELEMENT_TYPE_LINE) != 0 \
    && OGS_MAX_ELEMENT_DIM >= 0 && OGS_MAX_ELEMENT_ORDER >= 1
        _builder[std::type_index(typeid(MeshLib::Point))] =
                OGS_NEW_LOCAL_ASSEMBLER(NumLib::ShapePoint1);
#endif

#if (OGS_ENABLED_ELEMENTS & ENABLED_ELEMENT_TYPE_LINE) != 0 \
        && OGS_MAX_ELEMENT_DIM >= 1 && OGS_MAX_ELEMENT_ORDER >= 1
        _builder[std::type_index(typeid(MeshLib::Line))] =
            OGS_NEW_LOCAL_ASSEMBLER(NumLib::ShapeLine2);
#endif

#if (OGS_ENABLED_ELEMENTS & ENABLED_ELEMENT_TYPE_LINE) != 0 \
    && OGS_MAX_ELEMENT_DIM >= 1 && OGS_MAX_ELEMENT_ORDER >= 2
        _builder[std::type_index(typeid(MeshLib::Line3))] =
            OGS_NEW_LOCAL_ASSEMBLER(NumLib::ShapeLine3);
#endif


        // /// Quads and Hexahedra ///////////////////////////////////

#if (OGS_ENABLED_ELEMENTS & ENABLED_ELEMENT_TYPE_QUAD) != 0 \
    && OGS_MAX_ELEMENT_DIM >= 2 && OGS_MAX_ELEMENT_ORDER >= 1
        _builder[std::type_index(typeid(MeshLib::Quad))] =
            OGS_NEW_LOCAL_ASSEMBLER(NumLib::ShapeQuad4);
#endif

#if (OGS_ENABLED_ELEMENTS & ENABLED_ELEMENT_TYPE_CUBOID) != 0 \
    && OGS_MAX_ELEMENT_DIM >= 3 && OGS_MAX_ELEMENT_ORDER >= 1
        _builder[std::type_index(typeid(MeshLib::Hex))] =
            OGS_NEW_LOCAL_ASSEMBLER(NumLib::ShapeHex8);
#endif

#if (OGS_ENABLED_ELEMENTS & ENABLED_ELEMENT_TYPE_QUAD) != 0 \
    && OGS_MAX_ELEMENT_DIM >= 2 && OGS_MAX_ELEMENT_ORDER >= 2
        _builder[std::type_index(typeid(MeshLib::Quad8))] =
            OGS_NEW_LOCAL_ASSEMBLER(NumLib::ShapeQuad8);
        _builder[std::type_index(typeid(MeshLib::Quad9))] =
            OGS_NEW_LOCAL_ASSEMBLER(NumLib::ShapeQuad9);
#endif

#if (OGS_ENABLED_ELEMENTS & ENABLED_ELEMENT_TYPE_CUBOID) != 0 \
    && OGS_MAX_ELEMENT_DIM >= 3 && OGS_MAX_ELEMENT_ORDER >= 2
        _builder[std::type_index(typeid(MeshLib::Hex20))] =
            OGS_NEW_LOCAL_ASSEMBLER(NumLib::ShapeHex20);
#endif


        // /// Simplices ////////////////////////////////////////////////

#if (OGS_ENABLED_ELEMENTS & ENABLED_ELEMENT_TYPE_TRI) != 0 \
    && OGS_MAX_ELEMENT_DIM >= 2 && OGS_MAX_ELEMENT_ORDER >= 1
        _builder[std::type_index(typeid(MeshLib::Tri))] =
            OGS_NEW_LOCAL_ASSEMBLER(NumLib::ShapeTri3);
#endif

#if (OGS_ENABLED_ELEMENTS & ENABLED_ELEMENT_TYPE_SIMPLEX) != 0 \
    && OGS_MAX_ELEMENT_DIM >= 3 && OGS_MAX_ELEMENT_ORDER >= 1
        _builder[std::type_index(typeid(MeshLib::Tet))] =
            OGS_NEW_LOCAL_ASSEMBLER(NumLib::ShapeTet4);
#endif

#if (OGS_ENABLED_ELEMENTS & ENABLED_ELEMENT_TYPE_TRI) != 0 \
    && OGS_MAX_ELEMENT_DIM >= 2 && OGS_MAX_ELEMENT_ORDER >= 2
        _builder[std::type_index(typeid(MeshLib::Tri6))] =
            OGS_NEW_LOCAL_ASSEMBLER(NumLib::ShapeTri6);
#endif

#if (OGS_ENABLED_ELEMENTS & ENABLED_ELEMENT_TYPE_SIMPLEX) != 0 \
    && OGS_MAX_ELEMENT_DIM >= 3 && OGS_MAX_ELEMENT_ORDER >= 2
        _builder[std::type_index(typeid(MeshLib::Tet10))] =
            OGS_NEW_LOCAL_ASSEMBLER(NumLib::ShapeTet10);
#endif


        // /// Prisms ////////////////////////////////////////////////////

#if (OGS_ENABLED_ELEMENTS & ENABLED_ELEMENT_TYPE_PRISM) != 0 \
    && OGS_MAX_ELEMENT_DIM >= 3 && OGS_MAX_ELEMENT_ORDER >= 1
        _builder[std::type_index(typeid(MeshLib::Prism))] =
            OGS_NEW_LOCAL_ASSEMBLER(NumLib::ShapePrism6);
#endif

#if (OGS_ENABLED_ELEMENTS & ENABLED_ELEMENT_TYPE_PRISM) != 0 \
    && OGS_MAX_ELEMENT_DIM >= 3 && OGS_MAX_ELEMENT_ORDER >= 2
        _builder[std::type_index(typeid(MeshLib::Prism15))] =
            OGS_NEW_LOCAL_ASSEMBLER(NumLib::ShapePrism15);
#endif

        // /// Pyramids //////////////////////////////////////////////////

#if (OGS_ENABLED_ELEMENTS & ENABLED_ELEMENT_TYPE_PYRAMID) != 0 \
    && OGS_MAX_ELEMENT_DIM >= 3 && OGS_MAX_ELEMENT_ORDER >= 1
        _builder[std::type_index(typeid(MeshLib::Pyramid))] =
            OGS_NEW_LOCAL_ASSEMBLER(NumLib::ShapePyra5);
#endif

#if (OGS_ENABLED_ELEMENTS & ENABLED_ELEMENT_TYPE_PYRAMID) != 0 \
    && OGS_MAX_ELEMENT_DIM >= 3 && OGS_MAX_ELEMENT_ORDER >= 2
        _builder[std::type_index(typeid(MeshLib::Pyramid13))] =
            OGS_NEW_LOCAL_ASSEMBLER(NumLib::ShapePyra13);
#endif
    }

    /// Sets the provided data_ptr to the newly created local assembler data.
    void operator()(
            const MeshLib::Element& e,
            LocalAssemblerDataInterface_<GlobalMatrix_, GlobalVector_>*& data_ptr,
            std::size_t const local_matrix_size,
            unsigned const integration_order,
            ConstructorArgs&&... args)
    {
        auto const type_idx = std::type_index(typeid(e));
        auto it = _builder.find(type_idx);

        if (it != _builder.end()) {
            data_ptr = it->second(
                           e, local_matrix_size, integration_order,
                           std::forward<ConstructorArgs>(args)...);
        } else {
            ERR("You are trying to build a local assembler for an unknown mesh element type (%s)."
                " Maybe you have disabled this mesh element type in your build configuration.",
                type_idx.name());
            std::abort();
        }
    }

private:
    /// Mapping of element types to local assembler constructors.
    std::unordered_map<
        std::type_index,
        std::function<LocalAssemblerDataInterface_<GlobalMatrix_, GlobalVector_>*(
            MeshLib::Element const& e,
            std::size_t const local_matrix_size,
            unsigned const integration_order,
            ConstructorArgs&&...)>
        > _builder;
};

}   // namespace AssemblerLib


#undef ENABLED_ELEMENT_TYPE_SIMPLEX
#undef ENABLED_ELEMENT_TYPE_CUBOID
#undef ENABLED_ELEMENT_TYPE_PYRAMID
#undef ENABLED_ELEMENT_TYPE_PRISM
#undef ENABLED_ELEMENT_TYPE_LINE
#undef ENABLED_ELEMENT_TYPE_TRI
#undef ENABLED_ELEMENT_TYPE_QUAD
#undef OGS_ENABLED_ELEMENTS
#undef OGS_NEW_LOCAL_ASSEMBLER

#endif  // ASSEMBLER_LIB_LOCALDATAINITIALIZER_H_
