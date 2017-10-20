/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
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
#include "GlobalVectorFromNamedFunction.h"
#include "ProcessVariable.h"
#include "CoupledSolutionsForStaggeredScheme.h"

namespace ProcessLib
{
Process::Process(
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    unsigned const integration_order,
    std::vector<std::reference_wrapper<ProcessVariable>>&& process_variables,
    SecondaryVariableCollection&& secondary_variables,
    NumLib::NamedFunctionCaller&& named_function_caller)
    : _mesh(mesh),
      _secondary_variables(std::move(secondary_variables)),
      _named_function_caller(std::move(named_function_caller)),
      _global_assembler(std::move(jacobian_assembler)),
      _is_monolithic_scheme(true),
      _coupled_solutions(nullptr),
      _integration_order(integration_order),
      _process_variables(std::move(process_variables)),
      _boundary_conditions(parameters)
{
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
    _boundary_conditions.addBCsForProcessVariables(
        _process_variables, *_local_to_global_index_map, _integration_order);

    for (int variable_id = 0;
         variable_id < static_cast<int>(_process_variables.size());
         ++variable_id)
    {
        ProcessVariable& pv = _process_variables[variable_id];
        auto sts =
            pv.createSourceTerms(*_local_to_global_index_map, variable_id,
                                 _integration_order);

        std::move(sts.begin(), sts.end(),
                  std::back_inserter(_source_terms));
    }
}

void Process::setVariableInitialCondition(const int variable_id, double const t,
                                          GlobalVector& x)
{
    SpatialPosition pos;

    auto const& pv = _process_variables[variable_id];
    DBUG("Set the initial condition of variable %s.",
         pv.get().getName().data());

    auto const& ic = pv.get().getInitialCondition();

    auto const num_comp = pv.get().getNumberOfComponents();

    const int mesh_subset_id = _is_monolithic_scheme ? variable_id : 0;

    for (int component_id = 0; component_id < num_comp; ++component_id)
    {
        auto const& mesh_subsets = _local_to_global_index_map->getMeshSubsets(
            mesh_subset_id, component_id);
        for (auto const& mesh_subset : mesh_subsets)
        {
            auto const mesh_id = mesh_subset->getMeshID();
            for (auto const* node : mesh_subset->getNodes())
            {
                MeshLib::Location const l(mesh_id, MeshLib::MeshItemType::Node,
                                          node->getID());

                pos.setNodeID(node->getID());
                auto const& ic_value = ic(t, pos);

                auto global_index =
                    std::abs(_local_to_global_index_map->getGlobalIndex(
                        l, mesh_subset_id, component_id));
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

void Process::setInitialConditions(const unsigned pcs_id, double const t,
                                   GlobalVector& x)
{
    if (!_is_monolithic_scheme)
    {
        setVariableInitialCondition(pcs_id, t, x);
        return;
    }

    for (std::size_t variable_id = 0; variable_id < _process_variables.size();
         variable_id++)
    {
        setVariableInitialCondition(variable_id, t, x);
    }
}

MathLib::MatrixSpecifications Process::getMatrixSpecifications() const
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

    _boundary_conditions.applyNaturalBC(t, x, K, b);

    for (auto const& st : _source_terms)
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

    // TODO apply BCs to Jacobian.
    _boundary_conditions.applyNaturalBC(t, x, K, b);
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
    if (_is_monolithic_scheme)
    {
        // Collect the mesh subsets in a vector.
        for (ProcessVariable const& pv : _process_variables)
        {
            std::generate_n(
                std::back_inserter(all_mesh_subsets),
                pv.getNumberOfComponents(),
                [&]() {
                    return MeshLib::MeshSubsets{_mesh_subset_all_nodes.get()};
                });
        }

        // Create a vector of the number of variable components
        for (ProcessVariable const& pv : _process_variables)
            vec_var_n_components.push_back(pv.getNumberOfComponents());
    }
    else  // for staggered scheme
    {
        // Assuming that all equations of the coupled process use the same
        // element order. Other cases can be considered by overloading this
        // member function in the derived class.

        // Collect the mesh subsets in a vector.
        std::generate_n(
            std::back_inserter(all_mesh_subsets),
            _process_variables[0].get().getNumberOfComponents(),
            [&]() {
                return MeshLib::MeshSubsets{_mesh_subset_all_nodes.get()};
            });

        // Create a vector of the number of variable components.
        vec_var_n_components.push_back(
            _process_variables[0].get().getNumberOfComponents());
    }
    _local_to_global_index_map =
        std::make_unique<NumLib::LocalToGlobalIndexMap>(
            std::move(all_mesh_subsets), vec_var_n_components,
            NumLib::ComponentOrder::BY_LOCATION);
}

void Process::initializeExtrapolator()
{
    NumLib::LocalToGlobalIndexMap const* dof_table_single_component;
    bool manage_storage;

    if (_local_to_global_index_map->getNumberOfComponents() == 1 ||
        !_is_monolithic_scheme)
    {
        // For single-variable-single-component processes reuse the existing DOF
        // table.
        dof_table_single_component = _local_to_global_index_map.get();
        manage_storage = false;
    }
    else
    {
        // Otherwise construct a new DOF table.
        std::vector<MeshLib::MeshSubsets> all_mesh_subsets_single_component;
        all_mesh_subsets_single_component.emplace_back(
            _mesh_subset_all_nodes.get());

        dof_table_single_component = new NumLib::LocalToGlobalIndexMap(
            std::move(all_mesh_subsets_single_component),
            // by location order is needed for output
            NumLib::ComponentOrder::BY_LOCATION);
        manage_storage = true;
    }

    std::unique_ptr<NumLib::Extrapolator> extrapolator(
        new NumLib::LocalLinearLeastSquaresExtrapolator(
            *dof_table_single_component));

    // TODO Later on the DOF table can change during the simulation!
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

void Process::postTimestep(GlobalVector const& x)
{
    MathLib::LinAlg::setLocalAccessibleVector(x);
    postTimestepConcreteProcess(x);
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

void Process::setCoupledSolutionsOfPreviousTimeStep()
{
    const auto number_of_coupled_solutions =
        _coupled_solutions->coupled_xs.size();
    _coupled_solutions->coupled_xs_t0.reserve(number_of_coupled_solutions);
    for (std::size_t i = 0; i < number_of_coupled_solutions; i++)
    {
        const auto x_t0 = getPreviousTimeStepSolution(i);
        if (!x_t0)
        {
            OGS_FATAL(
                "Memory is not allocated for the global vector "
                "of the solution of the previous time step for the ."
                "staggered scheme.\n It can be done by overloading "
                "Process::preTimestepConcreteProcess"
                "(ref. HTProcess::preTimestepConcreteProcess) ");
        }

        MathLib::LinAlg::setLocalAccessibleVector(*x_t0);
        _coupled_solutions->coupled_xs_t0.emplace_back(x_t0);
    }
}

}  // namespace ProcessLib
