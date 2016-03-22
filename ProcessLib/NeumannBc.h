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
template <typename GlobalSetup>
class NeumannBc
{
public:
    using GlobalVector = typename GlobalSetup::VectorType;
    using GlobalMatrix = typename GlobalSetup::MatrixType;

    /// Create a Neumann boundary condition process from given config,
    /// DOF-table, and a mesh subset for a given variable and its component.
    /// A local DOF-table, a subset of the given one, is constructed.
    NeumannBc(
        NeumannBcConfig const& bc,
        unsigned const integration_order,
        AssemblerLib::LocalToGlobalIndexMap const& local_to_global_index_map,
        std::size_t const variable_id,
        std::size_t const component_id)
        : _function(*bc.getFunction()),
          _integration_order(integration_order)
    {
        assert(component_id < local_to_global_index_map.getNumComponents());

        // deep copy because the neumann bc config destroys the elements.
        std::transform(bc.elementsBegin(), bc.elementsEnd(),
                std::back_inserter(_elements),
                std::mem_fn(&MeshLib::Element::clone));

        std::vector<MeshLib::Node*> nodes = MeshLib::getUniqueNodes(_elements);

        auto const& mesh_subsets =
            local_to_global_index_map.getMeshSubsets(variable_id, component_id);

        // TODO extend the node intersection to all parts of mesh_subsets, i.e.
        // to each of the MeshSubset in the mesh_subsets.
        _mesh_subset_all_nodes =
            mesh_subsets.getMeshSubset(0).getIntersectionByNodes(nodes);
        std::unique_ptr<MeshLib::MeshSubsets> all_mesh_subsets{
            new MeshLib::MeshSubsets{_mesh_subset_all_nodes}};

        // Create local DOF table from intersected mesh subsets for the given
        // variable and component ids.
        _local_to_global_index_map.reset(
            local_to_global_index_map.deriveBoundaryConstrainedMap(
                variable_id, component_id, std::move(all_mesh_subsets),
                _elements));
    }

    ~NeumannBc()
    {
        delete _mesh_subset_all_nodes;

        for (auto e : _elements)
            delete e;

        for (auto p : _local_assemblers)
            delete p;
    }

    /// Calls local assemblers which calculate their contributions to the global
    /// matrix and the right-hand-side.
    void integrate(GlobalSetup const& global_setup,
                   const double t, GlobalVector& b)
    {
        global_setup.execute(*_global_assembler, _local_assemblers, t, b);
    }

    void initialize(GlobalSetup const& global_setup,
        unsigned global_dim)
    {
        if (global_dim==1)
            initialize<1u>(global_setup);
        else if (global_dim==2)
            initialize<2u>(global_setup);
        else if (global_dim==3)
            initialize<3u>(global_setup);
    }


private:
    /// Allocates the local assemblers for each element and stores references to
    /// global matrix and the right-hand-side.
    template <unsigned GlobalDim>
    void
    initialize(GlobalSetup const& global_setup)
    {
        // Shape matrices initializer
        using LocalDataInitializer = AssemblerLib::LocalDataInitializer<
            LocalNeumannBcAsmDataInterface,
            LocalNeumannBcAsmData,
            GlobalMatrix, GlobalVector,
            GlobalDim,
            std::function<double (MeshLib::Element const&)> const&>;

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
                _integration_order,
                elementValueLookup);

        DBUG("Create global assembler.");
        _global_assembler.reset(
            new GlobalAssembler(*_local_to_global_index_map));
    }


    /// The right-hand-side function of the Neumann boundary condition given as
    /// \f$ \alpha(x) \, \partial u(x) / \partial n = \text{_function}(x)\f$.
    MathLib::ConstantFunction<double> const _function;

    /// Vector of lower-dimensional elements on which the boundary condition is
    /// defined.
    std::vector<MeshLib::Element*> _elements;

    MeshLib::MeshSubset const* _mesh_subset_all_nodes = nullptr;

    /// Local dof table, a subset of the global one restricted to the
    /// participating #_elements of the boundary condition.
    std::unique_ptr<AssemblerLib::LocalToGlobalIndexMap> _local_to_global_index_map;

    /// Integration order for integration over the lower-dimensional elements of
    /// the #_function.
    unsigned const _integration_order;

    using GlobalAssembler =
        AssemblerLib::VectorMatrixAssembler<
            GlobalMatrix, GlobalVector,
            NumLib::ODESystemTag::NeumannBC>;

    std::unique_ptr<GlobalAssembler> _global_assembler;

    using LocalAssembler = LocalNeumannBcAsmDataInterface<
        GlobalMatrix, GlobalVector>;

    /// Local assemblers for each element of #_elements.
    std::vector<LocalAssembler*> _local_assemblers;

};

}   // namespace ProcessLib

#endif  // PROCESS_LIB_NEUMANN_BC_H_
