/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "TESProcess.h"

// TODO Copied from VectorMatrixAssembler. Could be provided by the DOF table.
inline NumLib::LocalToGlobalIndexMap::RowColumnIndices
getRowColumnIndices_(std::size_t const id,
                     NumLib::LocalToGlobalIndexMap const& dof_table,
                     std::vector<GlobalIndexType>& indices)
{
    assert(dof_table.size() > id);
    indices.clear();

    // Local matrices and vectors will always be ordered by component,
    // no matter what the order of the global matrix is.
    for (unsigned c = 0; c < dof_table.getNumComponents(); ++c)
    {
        auto const& idcs = dof_table(id, c).rows;
        indices.reserve(indices.size() + idcs.size());
        indices.insert(indices.end(), idcs.begin(), idcs.end());
    }

    return NumLib::LocalToGlobalIndexMap::RowColumnIndices(indices,
                                                                 indices);
}

template <typename GlobalVector>
void getVectorValues(
    GlobalVector const& x,
    NumLib::LocalToGlobalIndexMap::RowColumnIndices const& r_c_indices,
    std::vector<double>& local_x)
{
    auto const& indices = r_c_indices.rows;
    local_x.clear();
    local_x.reserve(indices.size());

    for (auto i : indices)
    {
        local_x.emplace_back(x.get(i));
    }
}

// TODO that essentially duplicates code which is also present in ProcessOutput.
template <typename GlobalVector>
double getNodalValue(GlobalVector const& x, MeshLib::Mesh const& mesh,
                     NumLib::LocalToGlobalIndexMap const& dof_table,
                     std::size_t const node_id,
                     std::size_t const global_component_id)
{
    MeshLib::Location const l{mesh.getID(), MeshLib::MeshItemType::Node,
                              node_id};

    auto const index = dof_table.getLocalIndex(
        l, global_component_id, x.getRangeBegin(), x.getRangeEnd());

    return x.get(index);
}

namespace ProcessLib
{
namespace TES
{

TESProcess::TESProcess(
    MeshLib::Mesh& mesh,
    typename Process::NonlinearSolver& nonlinear_solver,
    std::unique_ptr<typename Process::TimeDiscretization>&&
        time_discretization,
    std::vector<std::reference_wrapper<ProcessVariable>>&& process_variables,
    SecondaryVariableCollection&& secondary_variables,
    ProcessOutput&& process_output,
    const BaseLib::ConfigTree& config)
    : Process(
          mesh, nonlinear_solver, std::move(time_discretization),
          std::move(process_variables), std::move(secondary_variables),
          std::move(process_output))
{
    DBUG("Create TESProcess.");

    // physical parameters for local assembly
    {
        std::vector<std::pair<std::string, double*>> params{
            {"fluid_specific_heat_source",
             &_assembly_params.fluid_specific_heat_source},
            {"fluid_specific_isobaric_heat_capacity", &_assembly_params.cpG},
            {"solid_specific_heat_source",
             &_assembly_params.solid_specific_heat_source},
            {"solid_heat_conductivity", &_assembly_params.solid_heat_cond},
            {"solid_specific_isobaric_heat_capacity", &_assembly_params.cpS},
            {"tortuosity", &_assembly_params.tortuosity},
            {"diffusion_coefficient",
             &_assembly_params.diffusion_coefficient_component},
            {"porosity", &_assembly_params.poro},
            {"solid_density_dry", &_assembly_params.rho_SR_dry},
            {"solid_density_initial", &_assembly_params.initial_solid_density}};

        for (auto const& p : params)
        {
            if (auto const par = config.getConfigParameterOptional<double>(p.first))
            {
                DBUG("setting parameter `%s' to value `%g'", p.first.c_str(),
                     *par);
                *p.second = *par;
            }
        }
    }

    // characteristic values of primary variables
    {
        std::vector<std::pair<std::string, Trafo*>> const params{
            {"characteristic_pressure", &_assembly_params.trafo_p},
            {"characteristic_temperature", &_assembly_params.trafo_T},
            {"characteristic_vapour_mass_fraction", &_assembly_params.trafo_x}};

        for (auto const& p : params)
        {
            if (auto const par = config.getConfigParameterOptional<double>(p.first))
            {
                INFO("setting parameter `%s' to value `%g'", p.first.c_str(),
                     *par);
                *p.second = Trafo{*par};
            }
        }
    }

    // permeability
    if (auto par =
            config.getConfigParameterOptional<double>("solid_hydraulic_permeability"))
    {
        DBUG(
            "setting parameter `solid_hydraulic_permeability' to isotropic "
            "value `%g'",
            *par);
        const auto dim = mesh.getDimension();
        _assembly_params.solid_perm_tensor =
            Eigen::MatrixXd::Identity(dim, dim) * (*par);
    }

    // reactive system
    _assembly_params.react_sys = Adsorption::AdsorptionReaction::newInstance(
        config.getConfigSubtree("reactive_system"));

    // debug output
    if (auto const param =
            config.getConfigParameterOptional<bool>("output_element_matrices"))
    {
        DBUG("output_element_matrices: %s", (*param) ? "true" : "false");

        _assembly_params.output_element_matrices = *param;
    }

    // TODO somewhere else
    /*
    if (auto const param =
    config.getConfigParameterOptional<bool>("output_global_matrix"))
    {
        DBUG("output_global_matrix: %s", (*param) ? "true" : "false");

        this->_process_output.output_global_matrix = *param;
    }
    */
}


void TESProcess::initializeConcreteProcess(
    NumLib::LocalToGlobalIndexMap const& dof_table,
    MeshLib::Mesh const& mesh, unsigned const integration_order)
{
    DBUG("Create global assembler.");
    _global_assembler.reset(new GlobalAssembler(dof_table));

    ProcessLib::createLocalAssemblers<TESLocalAssembler>(
        mesh.getDimension(), mesh.getElements(), dof_table, integration_order,
        _local_assemblers, _assembly_params);

    // TODO move the two data members somewhere else.
    // for extrapolation of secondary variables
    std::vector<std::unique_ptr<MeshLib::MeshSubsets>>
        all_mesh_subsets_single_component;
    all_mesh_subsets_single_component.emplace_back(
        new MeshLib::MeshSubsets(this->_mesh_subset_all_nodes.get()));
    _local_to_global_index_map_single_component.reset(
        new NumLib::LocalToGlobalIndexMap(
            std::move(all_mesh_subsets_single_component),
            // by location order is needed for output
            NumLib::ComponentOrder::BY_LOCATION));

    {
        auto const& l = *_local_to_global_index_map_single_component;
        _extrapolator.reset(new ExtrapolatorImplementation(
            MathLib::MatrixSpecifications(l.dofSizeWithoutGhosts(),
                                          l.dofSizeWithoutGhosts(),
                                          &l.getGhostIndices(),
                                          nullptr),
            l));
    }

    // secondary variables
    auto add2nd = [&](std::string const& var_name, unsigned const n_components,
                      SecondaryVariableFunctions&& fcts) {
        this->_secondary_variables.addSecondaryVariable(var_name, n_components,
                                                        std::move(fcts));
    };
    auto makeEx =
        [&](TESIntPtVariables var) -> SecondaryVariableFunctions {
        return ProcessLib::makeExtrapolator(var, *_extrapolator,
                                            _local_assemblers);
    };

    add2nd("solid_density", 1, makeEx(TESIntPtVariables::SOLID_DENSITY));
    add2nd("reaction_rate", 1, makeEx(TESIntPtVariables::REACTION_RATE));
    add2nd("velocity_x", 1, makeEx(TESIntPtVariables::VELOCITY_X));
    if (mesh.getDimension() >= 2)
        add2nd("velocity_y", 1, makeEx(TESIntPtVariables::VELOCITY_Y));
    if (mesh.getDimension() >= 3)
        add2nd("velocity_z", 1, makeEx(TESIntPtVariables::VELOCITY_Z));

    add2nd("loading", 1, makeEx(TESIntPtVariables::LOADING));
    add2nd("reaction_damping_factor", 1,
           makeEx(TESIntPtVariables::REACTION_DAMPING_FACTOR));

    namespace PH = std::placeholders;
    using Self = TESProcess;

    add2nd("vapour_partial_pressure", 1,
           {std::bind(&Self::computeVapourPartialPressure, this, PH::_1, PH::_2,
                      PH::_3),
            nullptr});
    add2nd("relative_humidity", 1, {std::bind(&Self::computeRelativeHumidity,
                                              this, PH::_1, PH::_2, PH::_3),
                                    nullptr});
    add2nd("equilibrium_loading", 1,
           {std::bind(&Self::computeEquilibriumLoading, this, PH::_1, PH::_2,
                      PH::_3),
            nullptr});
}


void TESProcess::assembleConcreteProcess(const double t,
                                                      GlobalVector const& x,
                                                      GlobalMatrix& M,
                                                      GlobalMatrix& K,
                                                      GlobalVector& b)
{
    DBUG("Assemble TESProcess.");

    // Call global assembler for each local assembly item.
    ::detail::GlobalExecutorType::executeMemberDereferenced(*_global_assembler,
                                           &GlobalAssembler::assemble,
                                           _local_assemblers, t, x, M, K, b);
}


void TESProcess::preTimestep(GlobalVector const& x, const double t,
                                          const double delta_t)
{
    DBUG("new timestep");

    _assembly_params.delta_t = delta_t;
    _assembly_params.current_time = t;
    ++_assembly_params.timestep;  // TODO remove that

    _x_previous_timestep =
        MathLib::MatrixVectorTraits<GlobalVector>::newInstance(x);
}


void TESProcess::preIteration(const unsigned iter,
                                           GlobalVector const& /*x*/)
{
    _assembly_params.iteration_in_current_timestep = iter;
    ++_assembly_params.total_iteration;
    ++_assembly_params.number_of_try_of_iteration;
}


NumLib::IterationResult TESProcess::postIteration(
    GlobalVector const& x)
{
    if (this->_process_output.output_iteration_results)
    {
        DBUG("output results of iteration %li",
             _assembly_params.total_iteration);
        std::string fn =
            "tes_iter_" + std::to_string(_assembly_params.total_iteration) +
            +"_ts_" + std::to_string(_assembly_params.timestep) + "_" +
            std::to_string(_assembly_params.iteration_in_current_timestep) +
            "_" + std::to_string(_assembly_params.number_of_try_of_iteration) +
            ".vtu";

        this->output(fn, 0, x);
    }

    bool check_passed = true;

    if (!Trafo::constrained)
    {
        // bounds checking only has to happen if the vapour mass fraction is
        // non-logarithmic.

        std::vector<GlobalIndexType> indices_cache;
        std::vector<double> local_x_cache;
        std::vector<double> local_x_prev_ts_cache;

        auto check_variable_bounds = [&](std::size_t id,
                                         LocalAssembler& loc_asm) {
            auto const r_c_indices = getRowColumnIndices_(
                id, *this->_local_to_global_index_map, indices_cache);
            getVectorValues(x, r_c_indices, local_x_cache);
            getVectorValues(*_x_previous_timestep, r_c_indices,
                            local_x_prev_ts_cache);

            if (!loc_asm.checkBounds(local_x_cache, local_x_prev_ts_cache))
                check_passed = false;
        };

        ::detail::GlobalExecutorType::executeDereferenced(check_variable_bounds,
                                         _local_assemblers);
    }

    if (!check_passed)
        return NumLib::IterationResult::REPEAT_ITERATION;

    // TODO remove
    DBUG("ts %lu iteration %lu (in current ts: %lu) try %u accepted",
         _assembly_params.timestep, _assembly_params.total_iteration,
         _assembly_params.iteration_in_current_timestep,
         _assembly_params.number_of_try_of_iteration);

    _assembly_params.number_of_try_of_iteration = 0;

    return NumLib::IterationResult::SUCCESS;
}


typename TESProcess::GlobalVector const&
TESProcess::computeVapourPartialPressure(
    typename TESProcess::GlobalVector const& x,
    NumLib::LocalToGlobalIndexMap const& dof_table,
    std::unique_ptr<GlobalVector>& result_cache)
{
    assert(&dof_table == this->_local_to_global_index_map.get());

    auto const& dof_table_single = *_local_to_global_index_map_single_component;
    result_cache = MathLib::MatrixVectorTraits<GlobalVector>::newInstance(
        {dof_table_single.dofSizeWithoutGhosts(),
         dof_table_single.dofSizeWithoutGhosts(),
         &dof_table_single.getGhostIndices(), nullptr});

    GlobalIndexType nnodes = this->_mesh.getNumberOfNodes();

    for (GlobalIndexType node_id = 0; node_id < nnodes; ++node_id)
    {
        auto const p = getNodalValue(x, this->_mesh, dof_table, node_id,
                                     COMPONENT_ID_PRESSURE);
        auto const x_mV = getNodalValue(x, this->_mesh, dof_table, node_id,
                                        COMPONENT_ID_MASS_FRACTION);

        auto const x_nV = Adsorption::AdsorptionReaction::getMolarFraction(
            x_mV, _assembly_params.M_react, _assembly_params.M_inert);

        // TODO Problems with PETSc? (local vs. global index)
        result_cache->set(node_id, p * x_nV);
    }

    return *result_cache;
}


typename TESProcess::GlobalVector const&
TESProcess::computeRelativeHumidity(
    typename TESProcess::GlobalVector const& x,
    NumLib::LocalToGlobalIndexMap const& dof_table,
    std::unique_ptr<GlobalVector>& result_cache)
{
    assert(&dof_table == this->_local_to_global_index_map.get());

    auto const& dof_table_single = *_local_to_global_index_map_single_component;
    result_cache = MathLib::MatrixVectorTraits<GlobalVector>::newInstance(
        {dof_table_single.dofSizeWithoutGhosts(),
         dof_table_single.dofSizeWithoutGhosts(),
         &dof_table_single.getGhostIndices(), nullptr});

    GlobalIndexType nnodes = this->_mesh.getNumberOfNodes();

    for (GlobalIndexType node_id = 0; node_id < nnodes; ++node_id)
    {
        auto const p = getNodalValue(x, this->_mesh, dof_table, node_id,
                                     COMPONENT_ID_PRESSURE);
        auto const T = getNodalValue(x, this->_mesh, dof_table, node_id,
                                     COMPONENT_ID_TEMPERATURE);
        auto const x_mV = getNodalValue(x, this->_mesh, dof_table, node_id,
                                        COMPONENT_ID_MASS_FRACTION);

        auto const x_nV = Adsorption::AdsorptionReaction::getMolarFraction(
            x_mV, _assembly_params.M_react, _assembly_params.M_inert);

        auto const p_S =
            Adsorption::AdsorptionReaction::getEquilibriumVapourPressure(T);

        // TODO Problems with PETSc? (local vs. global index)
        result_cache->set(node_id, p * x_nV / p_S);
    }

    return *result_cache;
}


typename TESProcess::GlobalVector const&
TESProcess::computeEquilibriumLoading(
    typename TESProcess::GlobalVector const& x,
    NumLib::LocalToGlobalIndexMap const& dof_table,
    std::unique_ptr<GlobalVector>& result_cache)
{
    assert(&dof_table == this->_local_to_global_index_map.get());

    auto const& dof_table_single = *_local_to_global_index_map_single_component;
    result_cache = MathLib::MatrixVectorTraits<GlobalVector>::newInstance(
        {dof_table_single.dofSizeWithoutGhosts(),
         dof_table_single.dofSizeWithoutGhosts(),
         &dof_table_single.getGhostIndices(), nullptr});

    GlobalIndexType nnodes = this->_mesh.getNumberOfNodes();

    for (GlobalIndexType node_id = 0; node_id < nnodes; ++node_id)
    {
        auto const p = getNodalValue(x, this->_mesh, dof_table, node_id,
                                     COMPONENT_ID_PRESSURE);
        auto const T = getNodalValue(x, this->_mesh, dof_table, node_id,
                                     COMPONENT_ID_TEMPERATURE);
        auto const x_mV = getNodalValue(x, this->_mesh, dof_table, node_id,
                                        COMPONENT_ID_MASS_FRACTION);

        auto const x_nV = Adsorption::AdsorptionReaction::getMolarFraction(
            x_mV, _assembly_params.M_react, _assembly_params.M_inert);

        auto const p_V = p * x_nV;
        auto const C_eq =
            (p_V <= 0.0) ? 0.0
                         : _assembly_params.react_sys->getEquilibriumLoading(
                               p_V, T, _assembly_params.M_react);

        // TODO Problems with PETSc? (local vs. global index)
        result_cache->set(node_id, C_eq);
    }

    return *result_cache;
}

std::unique_ptr<TESProcess> createTESProcess(
    MeshLib::Mesh& mesh,
    typename Process::NonlinearSolver& nonlinear_solver,
    std::unique_ptr<typename Process::TimeDiscretization>&&
        time_discretization,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterBase>> const& /*parameters*/,
    BaseLib::ConfigTree const& config)
{
    config.checkConfigParameter("type", "TES");

    DBUG("Create TESProcess.");

    auto process_variables = findProcessVariables(
        variables, config,
        {"fluid_pressure", "temperature", "vapour_mass_fraction"});

    SecondaryVariableCollection
        secondary_variables{
            config.getConfigSubtreeOptional("secondary_variables"),
            {"solid_density", "reaction_rate", "velocity_x", "velocity_y",
             "velocity_z", "loading", "reaction_damping_factor",
             "vapour_partial_pressure", "relative_humidity",
             "equilibrium_loading"}};

    ProcessOutput process_output{
        config.getConfigSubtree("output"), process_variables,
        secondary_variables};

    return std::unique_ptr<TESProcess>{new TESProcess{
        mesh, nonlinear_solver, std::move(time_discretization),
        std::move(process_variables), std::move(secondary_variables),
        std::move(process_output), config}};
}

}  // namespace TES

}  // namespace ProcessLib
