/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "NeumannBc.h"

#include "NumLib/Fem/FiniteElement/TemplateIsoparametric.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"

#include "MeshLib/MeshSearch/NodeSearch.h"

#include "Utils/CreateLocalAssemblers.h"

namespace ProcessLib
{

NeumannBc::NeumannBc(
    NeumannBcConfig const& bc,
    unsigned const integration_order,
    NumLib::LocalToGlobalIndexMap const& local_to_global_index_map,
    int const variable_id,
    int const component_id)
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

NeumannBc::~NeumannBc()
{
    delete _mesh_subset_all_nodes;

    for (auto e : _elements)
        delete e;
}

void NeumannBc::integrate(const double t, GlobalVector& b)
{
    GlobalExecutor::executeMemberDereferenced(
                *_global_assembler, &GlobalAssembler::assemble,
                _local_assemblers, t, b);
}

void NeumannBc::initialize(unsigned global_dim)
{
    DBUG("Create global assembler.");
    _global_assembler.reset(
        new GlobalAssembler(*_local_to_global_index_map));

    auto elementValueLookup = [this](MeshLib::Element const&)
    {
        return _function();
    };

    createLocalAssemblers<LocalNeumannBcAsmData>(
        global_dim, _elements,
        *_local_to_global_index_map, _integration_order,
        _local_assemblers,
        elementValueLookup
        );
}

}   // namespace ProcessLib
