/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESSLIB_LIE_HYDROMECHANICS_LOCALDATAINITIALIZER_H_
#define PROCESSLIB_LIE_HYDROMECHANICS_LOCALDATAINITIALIZER_H_

#include <functional>
#include <memory>
#include <typeindex>
#include <typeinfo>
#include <type_traits>
#include <unordered_map>

#include "MeshLib/Elements/Elements.h"
#include "NumLib/DOF/LocalToGlobalIndexMap.h"
#include "NumLib/Fem/FiniteElement/LowerDimShapeTable.h"
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
template <
    typename LocalAssemblerInterface,
    template <typename, typename, typename, unsigned> class LocalAssemblerDataMatrix,
    template <typename, typename, typename, unsigned> class LocalAssemblerDataMatrixNearFracture,
    template <typename, typename, typename, unsigned> class LocalAssemblerDataFracture,
    unsigned GlobalDim, typename... ConstructorArgs>
class LocalDataInitializer final
{
public:
    using LADataIntfPtr = std::unique_ptr<LocalAssemblerInterface>;

    LocalDataInitializer(NumLib::LocalToGlobalIndexMap const& dof_table,
                         const unsigned shapefunction_order)
        : _dof_table(dof_table)
    {
        if (shapefunction_order != 2)
            OGS_FATAL(
                "The given shape function order %d is not supported.\nOnly "
                "shape functions of order 2 are supported.",
                shapefunction_order);
        // /// Quads and Hexahedra ///////////////////////////////////

#if (OGS_ENABLED_ELEMENTS & ENABLED_ELEMENT_TYPE_QUAD) != 0 && \
    OGS_MAX_ELEMENT_DIM >= 2 && OGS_MAX_ELEMENT_ORDER >= 2
        _builder[std::type_index(typeid(MeshLib::Quad8))] =
            makeLocalAssemblerBuilder<NumLib::ShapeQuad8>();
        _builder[std::type_index(typeid(MeshLib::Quad9))] =
            makeLocalAssemblerBuilder<NumLib::ShapeQuad9>();
#endif

        // /// Simplices ////////////////////////////////////////////////

#if (OGS_ENABLED_ELEMENTS & ENABLED_ELEMENT_TYPE_TRI) != 0 && \
    OGS_MAX_ELEMENT_DIM >= 2 && OGS_MAX_ELEMENT_ORDER >= 2
        _builder[std::type_index(typeid(MeshLib::Tri6))] =
            makeLocalAssemblerBuilder<NumLib::ShapeTri6>();
#endif

        // /// Lines ///////////////////////////////////
#if OGS_MAX_ELEMENT_DIM >= 3 && OGS_MAX_ELEMENT_ORDER >= 2
        _builder[std::type_index(typeid(MeshLib::Line3))] =
            makeLocalAssemblerBuilder<NumLib::ShapeLine3>();
#endif
    }

    /// Sets the provided \c data_ptr to the newly created local assembler data.
    ///
    /// \attention
    /// The index \c id is not necessarily the mesh item's id. Especially when
    /// having multiple meshes it will differ from the latter.
    void operator()(
            std::size_t const id,
            MeshLib::Element const& mesh_item,
            LADataIntfPtr& data_ptr,
            ConstructorArgs&&... args) const
    {
        auto const type_idx = std::type_index(typeid(mesh_item));
        auto const it = _builder.find(type_idx);

        if (it == _builder.end())
            OGS_FATAL("You are trying to build a local assembler for an unknown mesh element type (%s)."
                " Maybe you have disabled this mesh element type in your build configuration.",
                type_idx.name());

        auto const n_local_dof = _dof_table.getNumberOfElementDOF(id);
        auto const n_global_components = _dof_table.getNumberOfElementComponents(id);
        const std::vector<std::size_t> varIDs(_dof_table.getElementVariableIDs(id));

        std::vector<unsigned> dofIndex_to_localIndex;
        if (mesh_item.getDimension() == GlobalDim
            && n_global_components == GlobalDim + 1) // +1 for pressure
        {
            // normal matrix assembler
            dofIndex_to_localIndex.resize(n_local_dof);
            std::iota(dofIndex_to_localIndex.begin(), dofIndex_to_localIndex.end(), 0);
        }
        else
        {
            // matrix and fracture assemblers with enrichments
            dofIndex_to_localIndex.resize(n_local_dof);
            unsigned dof_id = 0;
            unsigned local_id = 0;
            std::vector<unsigned> vec_n_element_nodes;
            //TODO how to get the shape function order for each variable?
            vec_n_element_nodes.push_back(mesh_item.getNumberOfBaseNodes()); // pressure
            for (unsigned i=1; i<varIDs.size(); i++)
                vec_n_element_nodes.push_back(mesh_item.getNumberOfNodes()); // displacements

            for (unsigned i=0; i<varIDs.size(); i++)
            {
                auto const var_id = varIDs[i];
                auto const n_var_comp = _dof_table.getNumberOfVariableComponents(var_id);
                auto const n_var_element_nodes = vec_n_element_nodes[i];
                for (unsigned var_comp_id = 0; var_comp_id < n_var_comp; var_comp_id++)
                {
                    auto& mss = _dof_table.getMeshSubsets(var_id, var_comp_id);
                    assert(mss.size() == 1);
                    auto mesh_id = mss.getMeshSubset(0).getMeshID();
                    for (unsigned k=0; k<n_var_element_nodes; k++)
                    {
                        MeshLib::Location l(mesh_id,
                                            MeshLib::MeshItemType::Node,
                                            mesh_item.getNodeIndex(k));
                        auto global_index = _dof_table.getGlobalIndex(l, var_id, var_comp_id);
                        if (global_index != NumLib::MeshComponentMap::nop)
                            dofIndex_to_localIndex[dof_id++] = local_id;
                        local_id++;
                    }
                }
            }

        }

        data_ptr = it->second(
                       mesh_item, varIDs.size(), n_local_dof, dofIndex_to_localIndex,
                       std::forward<ConstructorArgs>(args)...);
    }

private:
    using LADataBuilder = std::function<LADataIntfPtr(
            MeshLib::Element const& e,
            std::size_t const n_variables,
            std::size_t const local_matrix_size,
            std::vector<unsigned> const& dofIndex_to_localIndex,
            ConstructorArgs&&...
        )>;

    template <typename ShapeFunctionDisplacement>
    using IntegrationMethod = typename NumLib::GaussIntegrationPolicy<
                typename ShapeFunctionDisplacement::MeshElement>::IntegrationMethod;

    template <typename ShapeFunctionDisplacement,
              typename ShapeFunctionPressure>
    using LADataMatrix =
        LocalAssemblerDataMatrix<ShapeFunctionDisplacement,
                                 ShapeFunctionPressure,
                                 IntegrationMethod<ShapeFunctionDisplacement>,
                                 GlobalDim>;

    template <typename ShapeFunctionDisplacement,
              typename ShapeFunctionPressure>
    using LADataMatrixNearFracture = LocalAssemblerDataMatrixNearFracture<
        ShapeFunctionDisplacement, ShapeFunctionPressure,
        IntegrationMethod<ShapeFunctionDisplacement>, GlobalDim>;

    template <typename ShapeFunctionDisplacement,
              typename ShapeFunctionPressure>
    using LAFractureData =
        LocalAssemblerDataFracture<ShapeFunctionDisplacement,
                                   ShapeFunctionPressure,
                                   IntegrationMethod<ShapeFunctionDisplacement>,
                                   GlobalDim>;

    /// A helper forwarding to the correct version of makeLocalAssemblerBuilder
    /// depending whether the global dimension is less than the shape function's
    /// dimension or not.
    template <typename ShapeFunctionDisplacement>
    static LADataBuilder makeLocalAssemblerBuilder()
    {
        return makeLocalAssemblerBuilder<ShapeFunctionDisplacement>(
            static_cast<std::integral_constant<
                bool, (GlobalDim >= ShapeFunctionDisplacement::DIM)>*>(
                nullptr));
    }

    /// Mapping of element types to local assembler constructors.
    std::unordered_map<std::type_index, LADataBuilder> _builder;

    NumLib::LocalToGlobalIndexMap const& _dof_table;

    // local assembler builder implementations.
private:
    /// Generates a function that creates a new LocalAssembler of type
    /// LAData<ShapeFunctionDisplacement>. Only functions with shape function's
    /// dimension less or equal to the global dimension are instantiated, e.g.
    /// following combinations of shape functions and global dimensions: (Line2,
    /// 1),
    /// (Line2, 2), (Line2, 3), (Hex20, 3) but not (Hex20, 2) or (Hex20, 1).
    template <typename ShapeFunctionDisplacement>
    static LADataBuilder makeLocalAssemblerBuilder(std::true_type*)
    {
        // (Lower order elements = Order(ShapeFunctionDisplacement) - 1).
        using ShapeFunctionPressure =
            typename NumLib::LowerDim<ShapeFunctionDisplacement>::type;
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
                        new LADataMatrix<ShapeFunctionDisplacement, ShapeFunctionPressure>{
                            e, n_variables, local_matrix_size, dofIndex_to_localIndex,
                            std::forward<ConstructorArgs>(args)...
                        }};
                } else {
                    return LADataIntfPtr{
                        new LADataMatrixNearFracture<ShapeFunctionDisplacement, ShapeFunctionPressure>{
                            e, n_variables, local_matrix_size, dofIndex_to_localIndex,
                            std::forward<ConstructorArgs>(args)...
                        }};
                }
            }
            return LADataIntfPtr{
                new LAFractureData<ShapeFunctionDisplacement, ShapeFunctionPressure>{
                    e, local_matrix_size, dofIndex_to_localIndex,
                    std::forward<ConstructorArgs>(args)...
                }};
        };
    }

    /// Returns nullptr for shape functions whose dimensions are less than the
    /// global dimension.
    template <typename ShapeFunctionDisplacement>
    static LADataBuilder makeLocalAssemblerBuilder(std::false_type*)
    {
        return nullptr;
    }
};

}  // namespace HydroMechanics
}  // namespace LIE
}  // namespace ProcessLib

#undef ENABLED_ELEMENT_TYPE_SIMPLEX
#undef ENABLED_ELEMENT_TYPE_CUBOID
#undef ENABLED_ELEMENT_TYPE_PYRAMID
#undef ENABLED_ELEMENT_TYPE_PRISM
#undef ENABLED_ELEMENT_TYPE_TRI
#undef ENABLED_ELEMENT_TYPE_QUAD
#undef OGS_ENABLED_ELEMENTS

#endif  // PROCESSLIB_LIE_HYDROMECHANICS_LOCALDATAINITIALIZER_H_
