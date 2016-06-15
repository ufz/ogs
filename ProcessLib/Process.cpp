/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "Process.h"

namespace ProcessLib
{
ProcessVariable& findProcessVariable(
    std::vector<ProcessVariable> const& variables,
    BaseLib::ConfigTree const& pv_config, std::string const& tag)
{
    // Find process variable name in process config.
    //! \ogs_file_special
    std::string const name = pv_config.getConfigParameter<std::string>(tag);

        // Find corresponding variable by name.
    auto variable = std::find_if(variables.cbegin(), variables.cend(),
                                 [&name](ProcessVariable const& v)
                                 {
                                     return v.getName() == name;
                                 });

    if (variable == variables.end())
    {
        OGS_FATAL(
            "Could not find process variable '%s' in the provided variables "
            "list for config tag <%s>.",
            name.c_str(), tag.c_str());
    }
    DBUG("Found process variable \'%s\' for config tag <%s>.",
         variable->getName().c_str(), tag.c_str());

    // Const cast is needed because of variables argument constness.
    return const_cast<ProcessVariable&>(*variable);
}

std::vector<std::reference_wrapper<ProcessVariable>>
findProcessVariables(
        std::vector<ProcessVariable> const& variables,
        BaseLib::ConfigTree const& process_config,
        std::initializer_list<std::string> tag_names)
{
    std::vector<std::reference_wrapper<ProcessVariable>> vars;
    vars.reserve(tag_names.size());

    //! \ogs_file_param{process__process_variables}
    auto const pv_conf = process_config.getConfigSubtree("process_variables");

    for (auto const& tag : tag_names) {
        vars.emplace_back(findProcessVariable(variables, pv_conf, tag));
    }

    return vars;
}

Process::Process(
    MeshLib::Mesh& mesh,
    NonlinearSolver& nonlinear_solver,
    std::unique_ptr<TimeDiscretization>&& time_discretization,
    std::vector<std::reference_wrapper<ProcessVariable>>&& process_variables,
    SecondaryVariableCollection&& secondary_variables,
    ProcessOutput&& process_output
    )
    : _mesh(mesh)
    , _secondary_variables(std::move(secondary_variables))
    , _process_output(std::move(process_output))
    , _nonlinear_solver(nonlinear_solver)
    , _time_discretization(std::move(time_discretization))
    , _process_variables(std::move(process_variables))
{}

void Process::initialize()
{
    DBUG("Initialize process.");

    DBUG("Construct dof mappings.");
    constructDofTable();

    DBUG("Compute sparsity pattern");
    computeSparsityPattern();

    initializeConcreteProcess(*_local_to_global_index_map, _mesh,
                              _integration_order);

    DBUG("Initialize boundary conditions.");

    // TODO That will only work with single component process variables!
    for (std::size_t global_component_id=0;
         global_component_id<_process_variables.size();
         ++global_component_id)
    {
        auto& pv = _process_variables[global_component_id];
        createDirichletBcs(pv, global_component_id);
        createNeumannBcs(pv, global_component_id);
    }

    for (auto& bc : _neumann_bcs)
        bc->initialize(_mesh.getDimension());
}

void Process::setInitialConditions(GlobalVector& x)
{
    DBUG("Set initial conditions.");

    // TODO That will only work with single component process variables!
    for (std::size_t global_component_id=0;
         global_component_id<_process_variables.size();
         ++global_component_id)
    {
        auto& pv = _process_variables[global_component_id];
        setInitialConditions(pv, global_component_id, x);
    }
}

void Process::constructDofTable()
{
    // Create single component dof in every of the mesh's nodes.
    _mesh_subset_all_nodes.reset(
        new MeshLib::MeshSubset(_mesh, &_mesh.getNodes()));

    // Collect the mesh subsets in a vector.
    std::vector<std::unique_ptr<MeshLib::MeshSubsets>> all_mesh_subsets;
    for (ProcessVariable const& pv : _process_variables)
    {
        std::generate_n(
            std::back_inserter(all_mesh_subsets),
            pv.getNumberOfComponents(),
            [&]()
            {
                return std::unique_ptr<MeshLib::MeshSubsets>{
                    new MeshLib::MeshSubsets{_mesh_subset_all_nodes.get()}};
            });
    }

    _local_to_global_index_map.reset(
        new NumLib::LocalToGlobalIndexMap(
            std::move(all_mesh_subsets),
            NumLib::ComponentOrder::BY_LOCATION));
}

void Process::setInitialConditions(ProcessVariable const& variable,
                          int const component_id,
                          GlobalVector& x)
{
    std::size_t const n_nodes = _mesh.getNumberOfNodes();
    for (std::size_t node_id = 0; node_id < n_nodes; ++node_id)
    {
        MeshLib::Location const l(_mesh.getID(),
                                  MeshLib::MeshItemType::Node, node_id);
        auto global_index = std::abs(
            _local_to_global_index_map->getGlobalIndex(l, component_id));
#ifdef USE_PETSC
        // The global indices of the ghost entries of the global
        // matrix or the global vectors need to be set as negative values
        // for equation assembly, however the global indices start from zero.
        // Therefore, any ghost entry with zero index is assigned an negative
        // value of the vector size or the matrix dimension.
        // To assign the initial value for the ghost entries,
        // the negative indices of the ghost entries are restored to zero.
        // checked hereby.
        if ( global_index == x.size() )
            global_index = 0;
#endif
        x.set(global_index,
              variable.getInitialConditionValue(*_mesh.getNode(node_id),
                                                component_id));
    }
}

void Process::createDirichletBcs(ProcessVariable& variable, int const component_id)
{
    MeshGeoToolsLib::MeshNodeSearcher& mesh_node_searcher =
        MeshGeoToolsLib::MeshNodeSearcher::getMeshNodeSearcher(
            variable.getMesh());

    variable.initializeDirichletBCs(std::back_inserter(_dirichlet_bcs),
                                    mesh_node_searcher,
                                    *_local_to_global_index_map,
                                    component_id);
}

void Process::createNeumannBcs(ProcessVariable& variable, int const component_id)
{
    // Find mesh nodes.
    MeshGeoToolsLib::MeshNodeSearcher& mesh_node_searcher =
        MeshGeoToolsLib::MeshNodeSearcher::getMeshNodeSearcher(
            variable.getMesh());
    MeshGeoToolsLib::BoundaryElementsSearcher mesh_element_searcher(
        variable.getMesh(), mesh_node_searcher);

    // Create a neumann BC for the process variable storing them in the
    // _neumann_bcs vector.
    variable.createNeumannBcs(
                              std::back_inserter(_neumann_bcs),
                              mesh_element_searcher,
                              _integration_order,
                              *_local_to_global_index_map,
                              0,  // 0 is the variable id TODO
                              component_id);
}

}  // namespace ProcessLib
