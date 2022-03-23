/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "Process.h"

#include "CoupledSolutionsForStaggeredScheme.h"
#include "NumLib/DOF/ComputeSparsityPattern.h"
#include "NumLib/Extrapolation/LocalLinearLeastSquaresExtrapolator.h"
#include "NumLib/ODESolver/ConvergenceCriterionPerComponent.h"
#include "ParameterLib/Parameter.h"
#include "ProcessVariable.h"

namespace ProcessLib
{
Process::Process(
    std::string name_,
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    unsigned const integration_order,
    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&&
        process_variables,
    SecondaryVariableCollection&& secondary_variables,
    const bool use_monolithic_scheme)
    : name(std::move(name_)),
      _mesh(mesh),
      _secondary_variables(std::move(secondary_variables)),
      _global_assembler(std::move(jacobian_assembler)),
      _use_monolithic_scheme(use_monolithic_scheme),
      _coupled_solutions(nullptr),
      _integration_order(integration_order),
      _process_variables(std::move(process_variables)),
      _boundary_conditions(
          [&](const std::size_t number_of_process_variables)
              -> std::vector<BoundaryConditionCollection>
          {
              std::vector<BoundaryConditionCollection> pcs_BCs;
              pcs_BCs.reserve(number_of_process_variables);
              for (std::size_t i = 0; i < number_of_process_variables; i++)
              {
                  pcs_BCs.emplace_back(BoundaryConditionCollection(parameters));
              }
              return pcs_BCs;
          }(_process_variables.size())),
      _source_term_collections(
          [&](const std::size_t number_of_processes)
              -> std::vector<SourceTermCollection>
          {
              std::vector<SourceTermCollection> pcs_sts;
              pcs_sts.reserve(number_of_processes);
              for (std::size_t i = 0; i < number_of_processes; i++)
              {
                  pcs_sts.emplace_back(SourceTermCollection(parameters));
              }
              return pcs_sts;
          }(_process_variables.size()))
{
}

void Process::initializeProcessBoundaryConditionsAndSourceTerms(
    const NumLib::LocalToGlobalIndexMap& dof_table, const int process_id)
{
    auto const& per_process_variables = _process_variables[process_id];
    auto& per_process_BCs = _boundary_conditions[process_id];

    per_process_BCs.addBCsForProcessVariables(per_process_variables, dof_table,
                                              _integration_order, *this);

    auto& per_process_sts = _source_term_collections[process_id];
    per_process_sts.addSourceTermsForProcessVariables(
        per_process_variables, dof_table, _integration_order);
}

void Process::initializeBoundaryConditions()
{
    // The number of processes is identical to the size of _process_variables,
    // the vector contains variables for different processes. See the
    // documentation of _process_variables.
    const std::size_t number_of_processes = _process_variables.size();
    for (std::size_t pcs_id = 0; pcs_id < number_of_processes; pcs_id++)
    {
        initializeProcessBoundaryConditionsAndSourceTerms(
            *_local_to_global_index_map, pcs_id);
    }
}

void Process::initialize()
{
    DBUG("Initialize process.");

    DBUG("Construct dof mappings.");
    constructDofTable();

    DBUG("Compute sparsity pattern");
    computeSparsityPattern();

    DBUG("Initialize the extrapolator");
    initializeExtrapolator();

    initializeConcreteProcess(*_local_to_global_index_map, _mesh,
                              _integration_order);

    DBUG("Initialize boundary conditions.");
    initializeBoundaryConditions();
}

void Process::setInitialConditions(
    std::vector<GlobalVector*>& process_solutions,
    std::vector<GlobalVector*> const& process_solutions_prev,
    double const t,
    int const process_id)
{
    auto& x = *process_solutions[process_id];
    auto& x_prev = *process_solutions_prev[process_id];

    // getDOFTableOfProcess can be overloaded by the specific process.
    auto const& dof_table_of_process = getDOFTable(process_id);

    auto const& per_process_variables = _process_variables[process_id];
    for (std::size_t variable_id = 0;
         variable_id < per_process_variables.size();
         variable_id++)
    {
        MathLib::LinAlg::setLocalAccessibleVector(x);
        ParameterLib::SpatialPosition pos;

        auto const& pv = per_process_variables[variable_id];
        DBUG("Set the initial condition of variable {:s} of process {:d}.",
             pv.get().getName().data(), process_id);

        auto const& ic = pv.get().getInitialCondition();

        auto const num_comp = pv.get().getNumberOfGlobalComponents();

        for (int component_id = 0; component_id < num_comp; ++component_id)
        {
            auto const& mesh_subset =
                dof_table_of_process.getMeshSubset(variable_id, component_id);
            auto const mesh_id = mesh_subset.getMeshID();
            for (auto const* node : mesh_subset.getNodes())
            {
                MeshLib::Location const l(mesh_id, MeshLib::MeshItemType::Node,
                                          node->getID());

                pos.setNodeID(node->getID());
                pos.setCoordinates(*node);
                auto const& ic_value = ic(t, pos);

                auto global_index =
                    std::abs(dof_table_of_process.getGlobalIndex(l, variable_id,
                                                                 component_id));
#ifdef USE_PETSC
                // The global indices of the ghost entries of the global
                // matrix or the global vectors need to be set as negative
                // values for equation assembly, however the global indices
                // start from zero. Therefore, any ghost entry with zero
                // index is assigned an negative value of the vector size
                // or the matrix dimension. To assign the initial value for
                // the ghost entries, the negative indices of the ghost
                // entries are restored to zero.
                if (global_index == x.size())
                    global_index = 0;
#endif
                x.set(global_index, ic_value[component_id]);
            }
        }
    }

    MathLib::LinAlg::finalizeAssembly(x);
    MathLib::LinAlg::copy(x, x_prev);  // pushState

    MathLib::LinAlg::setLocalAccessibleVector(x);
    MathLib::LinAlg::setLocalAccessibleVector(x_prev);

    setInitialConditionsConcreteProcess(process_solutions, t, process_id);
}

MathLib::MatrixSpecifications Process::getMatrixSpecifications(
    const int /*process_id*/) const
{
    auto const& l = *_local_to_global_index_map;
    return {l.dofSizeWithoutGhosts(), l.dofSizeWithoutGhosts(),
            &l.getGhostIndices(), &_sparsity_pattern};
}

void Process::updateDeactivatedSubdomains(double const time,
                                          const int process_id)
{
    auto const& variables_per_process = getProcessVariables(process_id);
    for (auto const& variable : variables_per_process)
    {
        variable.get().updateDeactivatedSubdomains(time);
    }
}

void Process::preAssemble(const double t, double const dt,
                          GlobalVector const& x)
{
    preAssembleConcreteProcess(t, dt, x);
}

void Process::assemble(const double t, double const dt,
                       std::vector<GlobalVector*> const& x,
                       std::vector<GlobalVector*> const& xdot,
                       int const process_id, GlobalMatrix& M, GlobalMatrix& K,
                       GlobalVector& b)
{
    assert(x.size() == xdot.size());
    for (std::size_t i = 0; i < x.size(); i++)
    {
        MathLib::LinAlg::setLocalAccessibleVector(*x[i]);
        MathLib::LinAlg::setLocalAccessibleVector(*xdot[i]);
    }

    assembleConcreteProcess(t, dt, x, xdot, process_id, M, K, b);

    // the last argument is for the jacobian, nullptr is for a unused jacobian
    _boundary_conditions[process_id].applyNaturalBC(t, x, process_id, K, b,
                                                    nullptr);

    // the last argument is for the jacobian, nullptr is for a unused jacobian
    _source_term_collections[process_id].integrate(t, *x[process_id], b,
                                                   nullptr);
}

void Process::assembleWithJacobian(const double t, double const dt,
                                   std::vector<GlobalVector*> const& x,
                                   std::vector<GlobalVector*> const& xdot,
                                   int const process_id, GlobalMatrix& M,
                                   GlobalMatrix& K, GlobalVector& b,
                                   GlobalMatrix& Jac)
{
    assert(x.size() == xdot.size());
    for (std::size_t i = 0; i < x.size(); i++)
    {
        MathLib::LinAlg::setLocalAccessibleVector(*x[i]);
        MathLib::LinAlg::setLocalAccessibleVector(*xdot[i]);
    }

    assembleWithJacobianConcreteProcess(t, dt, x, xdot, process_id, M, K, b,
                                        Jac);

    // TODO: apply BCs to Jacobian.
    _boundary_conditions[process_id].applyNaturalBC(t, x, process_id, K, b,
                                                    &Jac);

    // the last argument is for the jacobian, nullptr is for a unused jacobian
    _source_term_collections[process_id].integrate(t, *x[process_id], b, &Jac);
}

void Process::constructDofTable()
{
    if (_use_monolithic_scheme)
    {
        constructMonolithicProcessDofTable();

        return;
    }

    // For staggered scheme:
    const int specified_process_id = 0;
    constructDofTableOfSpecifiedProcessStaggeredScheme(specified_process_id);
}

void Process::constructMonolithicProcessDofTable()
{
    // Create single component dof in every of the mesh nodes.
    _mesh_subset_all_nodes =
        std::make_unique<MeshLib::MeshSubset>(_mesh, _mesh.getNodes());

    // Vector of mesh subsets.
    std::vector<MeshLib::MeshSubset> all_mesh_subsets;

    // Collect the mesh subsets in a vector for the components of each
    // variables.
    for (ProcessVariable const& pv : _process_variables[0])
    {
        std::generate_n(std::back_inserter(all_mesh_subsets),
                        pv.getNumberOfGlobalComponents(),
                        [&]() { return *_mesh_subset_all_nodes; });
    }

    // Create a vector of the number of variable components
    std::vector<int> vec_var_n_components;
    transform(cbegin(_process_variables[0]), cend(_process_variables[0]),
              back_inserter(vec_var_n_components),
              [](ProcessVariable const& pv)
              { return pv.getNumberOfGlobalComponents(); });

    _local_to_global_index_map =
        std::make_unique<NumLib::LocalToGlobalIndexMap>(
            std::move(all_mesh_subsets), vec_var_n_components,
            NumLib::ComponentOrder::BY_LOCATION);

    assert(_local_to_global_index_map);
}

void Process::constructDofTableOfSpecifiedProcessStaggeredScheme(
    const int specified_process_id)
{
    // Create single component dof in every of the mesh nodes.
    _mesh_subset_all_nodes =
        std::make_unique<MeshLib::MeshSubset>(_mesh, _mesh.getNodes());

    // Vector of mesh subsets.
    std::vector<MeshLib::MeshSubset> all_mesh_subsets;

    // Vector of the number of variable components
    std::vector<int> vec_var_n_components;
    // Collect the mesh subsets in a vector for each variables' components.
    std::generate_n(std::back_inserter(all_mesh_subsets),
                    _process_variables[specified_process_id][0]
                        .get()
                        .getNumberOfGlobalComponents(),
                    [&]() { return *_mesh_subset_all_nodes; });

    // Create a vector of the number of variable components.
    vec_var_n_components.push_back(_process_variables[specified_process_id][0]
                                       .get()
                                       .getNumberOfGlobalComponents());
    _local_to_global_index_map =
        std::make_unique<NumLib::LocalToGlobalIndexMap>(
            std::move(all_mesh_subsets), vec_var_n_components,
            NumLib::ComponentOrder::BY_LOCATION);

    assert(_local_to_global_index_map);
}

std::tuple<NumLib::LocalToGlobalIndexMap*, bool>
Process::getDOFTableForExtrapolatorData() const
{
    if (_local_to_global_index_map->getNumberOfGlobalComponents() == 1)
    {
        // For single-variable-single-component processes reuse the existing DOF
        // table.
        const bool manage_storage = false;
        return std::make_tuple(_local_to_global_index_map.get(),
                               manage_storage);
    }

    // Otherwise construct a new DOF table.
    std::vector<MeshLib::MeshSubset> all_mesh_subsets_single_component;
    all_mesh_subsets_single_component.emplace_back(*_mesh_subset_all_nodes);

    const bool manage_storage = true;

    return std::make_tuple(new NumLib::LocalToGlobalIndexMap(
                               std::move(all_mesh_subsets_single_component),
                               // by location order is needed for output
                               NumLib::ComponentOrder::BY_LOCATION),
                           manage_storage);
}

void Process::initializeExtrapolator()
{
    NumLib::LocalToGlobalIndexMap* dof_table_single_component;
    bool manage_storage;

    std::tie(dof_table_single_component, manage_storage) =
        getDOFTableForExtrapolatorData();

    std::unique_ptr<NumLib::Extrapolator> extrapolator(
        new NumLib::LocalLinearLeastSquaresExtrapolator(
            *dof_table_single_component));

    // TODO: Later on the DOF table can change during the simulation!
    _extrapolator_data = ExtrapolatorData(
        std::move(extrapolator), dof_table_single_component, manage_storage);
}

void Process::computeSparsityPattern()
{
    _sparsity_pattern =
        NumLib::computeSparsityPattern(*_local_to_global_index_map, _mesh);
}

void Process::preTimestep(std::vector<GlobalVector*> const& x, const double t,
                          const double delta_t, const int process_id)
{
    for (auto* const solution : x)
        MathLib::LinAlg::setLocalAccessibleVector(*solution);
    preTimestepConcreteProcess(x, t, delta_t, process_id);

    _boundary_conditions[process_id].preTimestep(t, x, process_id);
}

void Process::postTimestep(std::vector<GlobalVector*> const& x, const double t,
                           const double delta_t, int const process_id)
{
    for (auto* const solution : x)
        MathLib::LinAlg::setLocalAccessibleVector(*solution);
    postTimestepConcreteProcess(x, t, delta_t, process_id);

    _boundary_conditions[process_id].postTimestep(t, x, process_id);
}

void Process::postNonLinearSolver(GlobalVector const& x,
                                  GlobalVector const& xdot, const double t,
                                  double const dt, int const process_id)
{
    MathLib::LinAlg::setLocalAccessibleVector(x);
    MathLib::LinAlg::setLocalAccessibleVector(xdot);
    postNonLinearSolverConcreteProcess(x, xdot, t, dt, process_id);
}

void Process::computeSecondaryVariable(double const t,
                                       double const dt,
                                       std::vector<GlobalVector*> const& x,
                                       GlobalVector const& x_dot,
                                       int const process_id)
{
    for (auto const* solution : x)
        MathLib::LinAlg::setLocalAccessibleVector(*solution);
    MathLib::LinAlg::setLocalAccessibleVector(x_dot);

    computeSecondaryVariableConcrete(t, dt, x, x_dot, process_id);
}

void Process::preIteration(const unsigned iter, const GlobalVector& x)
{
    MathLib::LinAlg::setLocalAccessibleVector(x);
    preIterationConcreteProcess(iter, x);
}

NumLib::IterationResult Process::postIteration(const GlobalVector& x)
{
    MathLib::LinAlg::setLocalAccessibleVector(x);
    return postIterationConcreteProcess(x);
}

std::vector<GlobalIndexType>
Process::getIndicesOfResiduumWithoutInitialCompensation() const
{
    std::vector<GlobalIndexType> indices;

    for (std::size_t process_id = 0; process_id < _process_variables.size();
         process_id++)
    {
        auto const& dof_table_of_process = getDOFTable(process_id);

        auto const& per_process_variables = _process_variables[process_id];
        for (std::size_t variable_id = 0;
             variable_id < per_process_variables.size();
             variable_id++)
        {
            auto const& pv = per_process_variables[variable_id];
            DBUG("Set the initial condition of variable {:s} of process {:d}.",
                 pv.get().getName().data(), process_id);

            if ((pv.get().compensateNonEquilibriumInitialResiduum()))
            {
                continue;
            }

            auto const num_comp = pv.get().getNumberOfGlobalComponents();

            for (int component_id = 0; component_id < num_comp; ++component_id)
            {
                auto const& mesh_subset = dof_table_of_process.getMeshSubset(
                    variable_id, component_id);
                auto const mesh_id = mesh_subset.getMeshID();
                for (auto const* node : mesh_subset.getNodes())
                {
                    MeshLib::Location const l(
                        mesh_id, MeshLib::MeshItemType::Node, node->getID());

                    auto global_index =
                        std::abs(dof_table_of_process.getGlobalIndex(
                            l, variable_id, component_id));

                    indices.push_back(global_index);
                }
            }
        }
    }

    return indices;
}

}  // namespace ProcessLib
