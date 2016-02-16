/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESS_LIB_NEUMANN_BC_H_
#define PROCESS_LIB_NEUMANN_BC_H_

#include <memory>
#include <vector>

#include "NumLib/Fem/FiniteElement/TemplateIsoparametric.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"

#include "AssemblerLib/VectorMatrixAssembler.h"
#include "AssemblerLib/LocalToGlobalIndexMap.h"
#include "AssemblerLib/LocalDataInitializer.h"
#include "AssemblerLib/LocalAssemblerBuilder.h"
#include "MeshLib/MeshSubset.h"
#include "MeshLib/MeshSearch/NodeSearch.h"

#include "NeumannBcConfig.h"
#include "NeumannBcAssembler.h"

namespace ProcessLib
{

/// The NeumannBc class is a process which is integrating a single, Neumann type
/// boundary condition in to the global matrix and the right-hand-side.
///
/// The process operates on a set of elements and a subset of the DOF-table (the
/// local to global index map). For each element a local assembler is created
/// (NeumannBcAssembler).
///
/// In the construction phase the NeumannBcConfig together with global DOF-table
/// and mesh subset are used.
/// The creation of the local assemblers and binding to the global matrix and
/// right-hand-sides happen in the initialize() function.
/// The integration() function provides calls then the actual integration of the
/// Neumann boundary condition.
template <typename GlobalSetup_>
class NeumannBc
{
public:
    /// Create a Neumann boundary condition process from given config,
    /// DOF-table, and a mesh subset.
    /// A local DOF-table, a subset of the given one, is constructed.
    NeumannBc(
        NeumannBcConfig const& bc,
        unsigned const integration_order,
        AssemblerLib::LocalToGlobalIndexMap const& local_to_global_index_map,
        std::size_t const component_index,
        MeshLib::MeshSubset const& mesh_subset_all_nodes
        )
        :
          _function(*bc.getFunction()),
          _all_mesh_subsets(local_to_global_index_map.getNumComponents(), nullptr),
          _integration_order(integration_order)
    {
        assert(component_index < local_to_global_index_map.getNumComponents());

        // deep copy because the neumann bc config destroys the elements.
        std::transform(bc.elementsBegin(), bc.elementsEnd(),
                std::back_inserter(_elements),
                std::mem_fn(&MeshLib::Element::clone));

        std::vector<MeshLib::Node*> nodes = MeshLib::getUniqueNodes(_elements);

        _mesh_subset_all_nodes =
            mesh_subset_all_nodes.getIntersectionByNodes(nodes);

        _all_mesh_subsets[component_index] = new MeshLib::MeshSubsets(_mesh_subset_all_nodes);

        _local_to_global_index_map.reset(
            local_to_global_index_map.deriveBoundaryConstrainedMap(
                _all_mesh_subsets, _elements));
    }

    ~NeumannBc()
    {
        for (auto p : _all_mesh_subsets)
            delete p;

        delete _mesh_subset_all_nodes;

        for (auto e : _elements)
            delete e;

        for (auto p : _local_assemblers)
            delete p;
    }

    template <typename GlobalSetup>
    void
    initialize(GlobalSetup const& global_setup,
        unsigned global_dim)
    {
        if (global_dim==1)
            initialize<GlobalSetup, 1u>(global_setup);
        else if (global_dim==2)
            initialize<GlobalSetup, 2u>(global_setup);
        else if (global_dim==3)
            initialize<GlobalSetup, 3u>(global_setup);
    }

    /// Allocates the local assemblers for each element and stores references to
    /// global matrix and the right-hand-side.
    template <typename GlobalSetup, unsigned GlobalDim>
    void
    initialize(GlobalSetup const& global_setup)
    {
        // Shape matrices initializer
        using LocalDataInitializer = AssemblerLib::LocalDataInitializer<
            LocalNeumannBcAsmDataInterface,
            LocalNeumannBcAsmData,
            typename GlobalSetup::MatrixType,
            typename GlobalSetup::VectorType,
            GlobalDim>;

        LocalDataInitializer initializer;

        using LocalAssemblerBuilder =
            AssemblerLib::LocalAssemblerBuilder<
                MeshLib::Element,
                LocalDataInitializer>;

        // Populate the vector of local assemblers.
        _local_assemblers.resize(_elements.size());
        LocalAssemblerBuilder local_asm_builder(
            initializer, *_local_to_global_index_map);

        auto elementValueLookup = [this](MeshLib::Element const&)
        {
            return _function();
        };

        DBUG("Calling local Neumann assembler builder for Neumann boundary elements.");
        global_setup.transform(
                local_asm_builder,
                _elements,
                _local_assemblers,
                elementValueLookup,
                _integration_order);

        // TODO check
        DBUG("Create global assembler.");
        _global_assembler.reset(
            new GlobalAssembler(*_local_to_global_index_map));
    }

    /// Calls local assemblers which calculate their contributions to the global
    /// matrix and the right-hand-side.
    template <typename GlobalSetup>
    void
    integrate(GlobalSetup const& global_setup,
              typename GlobalSetup::VectorType const * x_curr = nullptr,
              typename GlobalSetup::VectorType const * x_prev_ts = nullptr
              )
    {
        _global_assembler->setX(x_curr, x_prev_ts);
        global_setup.execute(*_global_assembler, _local_assemblers);
    }


private:
    /// The right-hand-side function of the Neumann boundary condition given as
    /// \f$ \alpha(x) \, \partial u(x) / \partial n = \text{_function}(x)\f$.
    MathLib::ConstantFunction<double> const _function;

    /// Vector of lower-dimensional elements on which the boundary condition is
    /// defined.
    std::vector<MeshLib::Element*> _elements;

    MeshLib::MeshSubset const* _mesh_subset_all_nodes = nullptr;
    std::vector<MeshLib::MeshSubsets*> _all_mesh_subsets;

    /// Local dof table, a subset of the global one restricted to the
    /// participating #_elements of the boundary condition.
    std::unique_ptr<AssemblerLib::LocalToGlobalIndexMap> _local_to_global_index_map;

    /// Integration order for integration over the lower-dimensional elements of
    /// the #_function.
    unsigned const _integration_order;

    // TODO check
    using GlobalAssembler =
        AssemblerLib::VectorMatrixAssembler<
            typename GlobalSetup_::MatrixType,
            typename GlobalSetup_::VectorType,
            NumLib::ODESystemTag::FirstOrderImplicitQuasilinear>;

    std::unique_ptr<GlobalAssembler> _global_assembler;

    using LocalAssembler = LocalNeumannBcAsmDataInterface<
        typename GlobalSetup_::MatrixType, typename GlobalSetup_::VectorType>;

    /// Local assemblers for each element of #_elements.
    std::vector<LocalAssembler*> _local_assemblers;

};

}   // namespace ProcessLib

#endif  // PROCESS_LIB_NEUMANN_BC_H_
