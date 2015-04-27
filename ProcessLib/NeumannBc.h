/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
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

#include "NeumannBcConfig.h"
#include "NeumannBcAssembler.h"
#include "MeshLib/MeshSearcher.h"

namespace ProcessLib
{

template <typename GlobalSetup_>
class NeumannBc
{
public:
    NeumannBc(
        NeumannBcConfig* bc,
        unsigned const integration_order,
        AssemblerLib::LocalToGlobalIndexMap const& local_to_global_index_map,
        MeshLib::MeshSubset const& mesh_subset_all_nodes
        )
        :
          _function(*bc->getFunction()),
          _integration_order(integration_order)
    {
        // deep copy because the neumann bc config destroys the elements.
        std::transform(bc->elementsBegin(), bc->elementsEnd(),
                std::back_inserter(_elements),
                std::mem_fn(&MeshLib::Element::clone));

        std::vector<MeshLib::Node*> nodes = MeshLib::selectNodes(_elements);

        _mesh_subset_all_nodes =
            mesh_subset_all_nodes.getIntersectionByNodes(nodes);
        _all_mesh_subsets.push_back(new MeshLib::MeshSubsets(_mesh_subset_all_nodes));

        _local_to_global_index_map.reset(
            local_to_global_index_map.deriveBoundaryConstrainedMap(
                _all_mesh_subsets, _elements));
    }

    ~NeumannBc()
    {
        for (auto e : _elements)
            delete e;

        for (auto p : _local_assemblers)
            delete p;
    }

    template <typename GlobalSetup>
    void
    initialize(GlobalSetup const& global_setup,
        typename GlobalSetup::MatrixType& A,
        typename GlobalSetup::VectorType& rhs)
    {
        // Shape matrices initializer
        using LocalDataInitializer = AssemblerLib::LocalDataInitializer<
            LocalNeumannBcAsmDataInterface,
            LocalNeumannBcAsmData,
            typename GlobalSetup::MatrixType,
            typename GlobalSetup::VectorType>;

        LocalDataInitializer initializer;

        using LocalAssemblerBuilder =
            AssemblerLib::LocalAssemblerBuilder<
                MeshLib::Element,
                LocalDataInitializer>;

        // Populate the vector of local assemblers.
        _local_assemblers.resize(_elements.size());
        LocalAssemblerBuilder local_asm_builder(
            initializer, *_local_to_global_index_map);

        auto elementValueLookup = [this](MeshLib::Element const& e)
        {
            return _function();
        };

        DBUG("Calling local Neumann assembler builder for neumann boundary elements.");
        global_setup.execute(
                local_asm_builder,
                _elements,
                _local_assemblers,
                elementValueLookup,
                _integration_order);

        DBUG("Create global assembler.");
        _global_assembler.reset(
            new GlobalAssembler(A, rhs, *_local_to_global_index_map));
    }

    template <typename GlobalSetup>
    void
    integrate(GlobalSetup const& global_setup)
    {
        global_setup.execute(*_global_assembler, _local_assemblers);
    }


private:
    MathLib::ConstantFunction<double> const _function;

    std::vector<MeshLib::Element const*> _elements;
    MeshLib::MeshSubset const* _mesh_subset_all_nodes = nullptr;
    std::vector<MeshLib::MeshSubsets*> _all_mesh_subsets;

    std::unique_ptr<AssemblerLib::LocalToGlobalIndexMap> _local_to_global_index_map;

    unsigned const _integration_order;

    using GlobalAssembler =
        AssemblerLib::VectorMatrixAssembler<
            typename GlobalSetup_::MatrixType,
            typename GlobalSetup_::VectorType>;

    std::unique_ptr<GlobalAssembler> _global_assembler;

    using LocalAssembler = LocalNeumannBcAsmDataInterface<
        typename GlobalSetup_::MatrixType, typename GlobalSetup_::VectorType>;

    std::vector<LocalAssembler*> _local_assemblers;

};

}   // namespace ProcessLib

#endif  // PROCESS_LIB_NEUMANN_BC_H_
