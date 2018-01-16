/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "Process.h"

#include "BaseLib/Functional.h"
#include "NumLib/DOF/ComputeSparsityPattern.h"
#include "NumLib/Extrapolation/LocalLinearLeastSquaresExtrapolator.h"
#include "NumLib/ODESolver/ConvergenceCriterionPerComponent.h"
#include "ProcessLib/Output/GlobalVectorFromNamedFunction.h"

#include "ProcessVariable.h"
#include "CoupledSolutionsForStaggeredScheme.h"

namespace ProcessLib
{
Process::Process(
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    unsigned const integration_order,
    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&&
        process_variables,
    SecondaryVariableCollection&& secondary_variables,
    NumLib::NamedFunctionCaller&& named_function_caller,
    const bool use_monolithic_scheme)
    : _mesh(mesh),
      _secondary_variables(std::move(secondary_variables)),
      _named_function_caller(std::move(named_function_caller)),
      _global_assembler(std::move(jacobian_assembler)),
      _use_monolithic_scheme(use_monolithic_scheme),
      _coupled_solutions(nullptr),
      _integration_order(integration_order),
      _process_variables(std::move(process_variables)),
      _boundary_conditions([&](const std::size_t number_of_processes)
                               -> std::vector<BoundaryConditionCollection> {
          std::vector<BoundaryConditionCollection> pcs_BCs;
          pcs_BCs.reserve(number_of_processes);
          for (std::size_t i = 0; i < number_of_processes; i++)
          {
              pcs_BCs.emplace_back(BoundaryConditionCollection(parameters));
          }
          return pcs_BCs;
      }(_process_variables.size()))
{
}

void Process::initializeProcessBoundaryCondition(
    const NumLib::LocalToGlobalIndexMap& dof_table, const int process_id)
{
    auto const& per_process_variables = _process_variables[process_id];
    auto& per_process_BCs = _boundary_conditions[process_id];

    per_process_BCs.addBCsForProcessVariables(per_process_variables, dof_table,
                                              _integration_order);

    std::vector<std::unique_ptr<NodalSourceTerm>> per_process_source_terms;
    for (auto& pv : per_process_variables)
    {
        auto sts = pv.get().createSourceTerms(dof_table, 0, _integration_order);

        std::move(sts.begin(), sts.end(),
                  std::back_inserter(per_process_source_terms));
    }
    _source_terms.push_back(std::move(per_process_source_terms));
}

void Process::initializeBoundaryConditions()
{
    // The number of processes is identical to the size of _process_variables,
    // the vector contains variables for different processes. See the
    // documentation of _process_variables.
    const std::size_t number_of_processes = _process_variables.size();
    for (std::size_t pcs_id = 0; pcs_id < number_of_processes; pcs_id++)
    {
        initializeProcessBoundaryCondition(*_local_to_global_index_map, pcs_id);
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

    finishNamedFunctionsInitialization();

    DBUG("Initialize boundary conditions.");
    initializeBoundaryConditions();
}

void Process::setInitialConditions(const int process_id, double const t,
                                   GlobalVector& x)
{
    // getDOFTableOfProcess can be overloaded by the specific process.
    auto const dof_table_of_process = getDOFTable(process_id);

    auto const& per_process_variables = _process_variables[process_id];
    for (std::size_t variable_id = 0;
         variable_id < per_process_variables.size();
         variable_id++)
    {
        SpatialPosition pos;

        auto const& pv = per_process_variables[variable_id];
        DBUG("Set the initial condition of variable %s of process %d.",
             pv.get().getName().data(), process_id);

        auto const& ic = pv.get().getInitialCondition();

        auto const num_comp = pv.get().getNumberOfComponents();

        for (int component_id = 0; component_id < num_comp; ++component_id)
        {
            auto const& mesh_subsets =
                dof_table_of_process.getMeshSubsets(variable_id, component_id);
            for (auto const& mesh_subset : mesh_subsets)
            {
                auto const mesh_id = mesh_subset->getMeshID();
                for (auto const* node : mesh_subset->getNodes())
                {
                    MeshLib::Location const l(
                        mesh_id, MeshLib::MeshItemType::Node, node->getID());

                    pos.setNodeID(node->getID());
                    auto const& ic_value = ic(t, pos);

                    auto global_index =
                        std::abs(dof_table_of_process.getGlobalIndex(
                            l, variable_id, component_id));
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
    }
}

MathLib::MatrixSpecifications Process::getMatrixSpecifications(
    const int /*process_id*/) const
{
    auto const& l = *_local_to_global_index_map;
    return {l.dofSizeWithoutGhosts(), l.dofSizeWithoutGhosts(),
            &l.getGhostIndices(), &_sparsity_pattern};
}

void Process::preAssemble(const double t, GlobalVector const& x)
{
    preAssembleConcreteProcess(t, x);
}

void Process::assemble(const double t, GlobalVector const& x, GlobalMatrix& M,
                       GlobalMatrix& K, GlobalVector& b)
{
    MathLib::LinAlg::setLocalAccessibleVector(x);

    assembleConcreteProcess(t, x, M, K, b);

    const auto pcs_id =
        (_coupled_solutions) != nullptr ? _coupled_solutions->process_id : 0;
    _boundary_conditions[pcs_id].applyNaturalBC(t, x, K, b);

    auto& source_terms_per_pcs = _source_terms[pcs_id];
    for (auto& st : source_terms_per_pcs)
    {
        st->integrateNodalSourceTerm(t, b);
    }
}

void Process::assembleWithJacobian(const double t, GlobalVector const& x,
                                   GlobalVector const& xdot,
                                   const double dxdot_dx, const double dx_dx,
                                   GlobalMatrix& M, GlobalMatrix& K,
                                   GlobalVector& b, GlobalMatrix& Jac)
{
    MathLib::LinAlg::setLocalAccessibleVector(x);
    MathLib::LinAlg::setLocalAccessibleVector(xdot);

    assembleWithJacobianConcreteProcess(t, x, xdot, dxdot_dx, dx_dx, M, K, b,
                                        Jac);

    // TODO: apply BCs to Jacobian.
    const auto pcs_id =
        (_coupled_solutions) != nullptr ? _coupled_solutions->process_id : 0;
    _boundary_conditions[pcs_id].applyNaturalBC(t, x, K, b);
}

void Process::constructDofTable()
{
    // Create single component dof in every of the mesh's nodes.
    _mesh_subset_all_nodes =
        std::make_unique<MeshLib::MeshSubset>(_mesh, &_mesh.getNodes());

    // Vector of mesh subsets.
    std::vector<MeshLib::MeshSubsets> all_mesh_subsets;

    // Vector of the number of variable components
    std::vector<int> vec_var_n_components;
    if (_use_monolithic_scheme)
    {
        // Collect the mesh subsets in a vector.
        for (ProcessVariable const& pv : _process_variables[0])
        {
            std::generate_n(
                std::back_inserter(all_mesh_subsets),
                pv.getNumberOfComponents(),
                [&]() {
                    return MeshLib::MeshSubsets{_mesh_subset_all_nodes.get()};
                });
        }

        // Create a vector of the number of variable components
        for (ProcessVariable const& pv : _process_variables[0])
        {
            vec_var_n_components.push_back(pv.getNumberOfComponents());
        }
    }
    else  // for staggered scheme
    {
        // Assuming that all equations of the coupled process use the same
        // element order. Other cases can be considered by overloading this
        // member function in the derived class.

        // Collect the mesh subsets in a vector.
        std::generate_n(
            std::back_inserter(all_mesh_subsets),
            _process_variables[0][0].get().getNumberOfComponents(),
            [&]() {
                return MeshLib::MeshSubsets{_mesh_subset_all_nodes.get()};
            });

        // Create a vector of the number of variable components.
        vec_var_n_components.push_back(
            _process_variables[0][0].get().getNumberOfComponents());
    }
    _local_to_global_index_map =
        std::make_unique<NumLib::LocalToGlobalIndexMap>(
            std::move(all_mesh_subsets), vec_var_n_components,
            NumLib::ComponentOrder::BY_LOCATION);
}

std::tuple<NumLib::LocalToGlobalIndexMap*, bool>
Process::getDOFTableForExtrapolatorData() const
{
    if (_local_to_global_index_map->getNumberOfComponents() == 1)
    {
        // For single-variable-single-component processes reuse the existing DOF
        // table.
        const bool manage_storage = false;
        return std::make_tuple(_local_to_global_index_map.get(),
                               manage_storage);
    }

    // Otherwise construct a new DOF table.
    std::vector<MeshLib::MeshSubsets> all_mesh_subsets_single_component;
    all_mesh_subsets_single_component.emplace_back(
        _mesh_subset_all_nodes.get());

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

void Process::finishNamedFunctionsInitialization()
{
    _named_function_caller.applyPlugs();

    for (auto const& named_function :
         _named_function_caller.getNamedFunctions())
    {
        auto const& name = named_function.getName();
        _secondary_variables.addSecondaryVariable(
            name,
            {1, BaseLib::easyBind(
                    &GlobalVectorFromNamedFunction::call,
                    GlobalVectorFromNamedFunction(
                        _named_function_caller.getSpecificFunctionCaller(name),
                        _mesh, getSingleComponentDOFTable(),
                        _secondary_variable_context)),
             nullptr});
    }
}

void Process::computeSparsityPattern()
{
    _sparsity_pattern =
        NumLib::computeSparsityPattern(*_local_to_global_index_map, _mesh);
}

void Process::preTimestep(GlobalVector const& x, const double t,
                          const double delta_t, const int process_id)
{
    for (auto& cached_var : _cached_secondary_variables)
    {
        cached_var->setTime(t);
    }

    MathLib::LinAlg::setLocalAccessibleVector(x);
    preTimestepConcreteProcess(x, t, delta_t, process_id);
}

void Process::postTimestep(GlobalVector const& x, int const process_id)
{
    MathLib::LinAlg::setLocalAccessibleVector(x);
    postTimestepConcreteProcess(x, process_id);
}

void Process::postNonLinearSolver(GlobalVector const& x, const double t,
                                  int const process_id)
{
    MathLib::LinAlg::setLocalAccessibleVector(x);
    postNonLinearSolverProcess(x, t, process_id);
}

void Process::computeSecondaryVariable(const double t, GlobalVector const& x)
{
    MathLib::LinAlg::setLocalAccessibleVector(x);

    computeSecondaryVariableConcrete(t, x);
}

void Process::preIteration(const unsigned iter, const GlobalVector& x)
{
    // In every new iteration cached values of secondary variables are expired.
    for (auto& cached_var : _cached_secondary_variables)
    {
        cached_var->updateCurrentSolution(x, *_local_to_global_index_map);
    }

    MathLib::LinAlg::setLocalAccessibleVector(x);
    preIterationConcreteProcess(iter, x);
}

NumLib::IterationResult Process::postIteration(const GlobalVector& x)
{
    MathLib::LinAlg::setLocalAccessibleVector(x);
    return postIterationConcreteProcess(x);
}

}  // namespace ProcessLib
