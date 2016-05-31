/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESS_LIB_PROCESS_H_
#define PROCESS_LIB_PROCESS_H_

#include "MeshLib/IO/VtkIO/VtuInterface.h"
#include "NumLib/DOF/ComputeSparsityPattern.h"
#include "NumLib/ODESolver/ODESystem.h"
#include "NumLib/ODESolver/TimeDiscretization.h"
#include "NumLib/ODESolver/NonlinearSolver.h"

#include "DirichletBc.h"
#include "NeumannBc.h"
#include "NeumannBcAssembler.h"
#include "Parameter.h"
#include "ProcessOutput.h"
#include "ProcessVariable.h"
#include "UniformDirichletBoundaryCondition.h"

namespace MeshLib
{
class Mesh;
}

namespace ProcessLib
{

template <typename GlobalSetup>
class Process
        : public NumLib::ODESystem<typename GlobalSetup::MatrixType,
                                   typename GlobalSetup::VectorType,
                                   // TODO: later on use a simpler ODE system
                                   NumLib::ODESystemTag::FirstOrderImplicitQuasilinear,
                                   NumLib::NonlinearSolverTag::Newton>
{
public:
    using GlobalVector = typename GlobalSetup::VectorType;
    using GlobalMatrix = typename GlobalSetup::MatrixType;
    using Index = typename GlobalMatrix::IndexType;
    using NonlinearSolver = NumLib::NonlinearSolverBase<GlobalMatrix, GlobalVector>;
    using TimeDiscretization = NumLib::TimeDiscretization<GlobalVector>;

    Process(
        MeshLib::Mesh& mesh,
        NonlinearSolver& nonlinear_solver,
        std::unique_ptr<TimeDiscretization>&& time_discretization,
        std::vector<std::reference_wrapper<ProcessVariable>>&& process_variables,
        SecondaryVariableCollection<GlobalVector>&& secondary_variables,
        ProcessOutput<GlobalVector>&& process_output
        )
        : _mesh(mesh)
        , _secondary_variables(std::move(secondary_variables))
        , _process_output(std::move(process_output))
        , _nonlinear_solver(nonlinear_solver)
        , _time_discretization(std::move(time_discretization))
        , _process_variables(std::move(process_variables))
    {}

    /// Preprocessing before starting assembly for new timestep.
    virtual void preTimestep(GlobalVector const& /*x*/,
                             const double /*t*/, const double /*delta_t*/) {}

    /// Postprocessing after a complete timestep.
    virtual void postTimestep(GlobalVector const& /*x*/) {}

    /// Process output.
    /// The file_name is indicating the name of possible output file.
    void output(std::string const& file_name,
                const unsigned /*timestep*/,
                GlobalVector const& x) const
    {
        doProcessOutput(file_name, x, _mesh, *_local_to_global_index_map,
                        _process_variables, _secondary_variables, _process_output);
    }

    void initialize()
    {
        DBUG("Initialize process.");

        DBUG("Construct dof mappings.");
        constructDofTable();

        DBUG("Compute sparsity pattern");
        computeSparsityPattern();

        initializeConcreteProcess(*_local_to_global_index_map, _mesh,
                                  _integration_order);

        DBUG("Initialize boundary conditions.");
        for (int variable_id = 0;
             variable_id < static_cast<int>(_process_variables.size());
             ++variable_id)
        {
            ProcessVariable& pv = _process_variables[variable_id];
            for (int component_id = 0;
                 component_id < pv.getNumberOfComponents();
                 ++component_id)
            {
                createDirichletBcs(pv, variable_id, component_id);
                createNeumannBcs(pv, variable_id, component_id);
            }
        }

        for (auto& bc : _neumann_bcs)
            bc->initialize(_mesh.getDimension());
    }

    void setInitialConditions(GlobalVector& x)
    {
        DBUG("Set initial conditions.");
        for (int variable_id = 0;
             variable_id < static_cast<int>(_process_variables.size());
             ++variable_id)
        {
            ProcessVariable& pv = _process_variables[variable_id];
            for (int component_id = 0;
                 component_id < pv.getNumberOfComponents();
                 ++component_id)
            {
                setInitialConditions(pv, variable_id, component_id, x);
            }
        }
    }

    MathLib::MatrixSpecifications getMatrixSpecifications() const override final
    {
        auto const& l = *_local_to_global_index_map;
        return {l.dofSizeWithoutGhosts(), l.dofSizeWithoutGhosts(),
                &l.getGhostIndices(), &_sparsity_pattern};
    }

    void assemble(const double t, GlobalVector const& x,
                  GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b) override final
    {
        assembleConcreteProcess(t, x, M, K, b);

        // Call global assembler for each Neumann boundary local assembler.
        for (auto const& bc : _neumann_bcs)
            bc->integrate(t, b);
    }

    void assembleJacobian(
        const double t, GlobalVector const& x, GlobalVector const& xdot,
        const double dxdot_dx, GlobalMatrix const& M,
        const double dx_dx, GlobalMatrix const& K,
        GlobalMatrix& Jac) override final
    {
        assembleJacobianConcreteProcess(t, x, xdot, dxdot_dx, M, dx_dx, K, Jac);

        // TODO In this method one could check if the user wants to use an
        //      analytical or a numerical Jacobian. Then the right
        //      assembleJacobianConcreteProcess() method will be chosen.
        //      Additionally in the default implementation of said method one
        //      could provide a fallback to a numerical Jacobian. However, that
        //      would be in a sense implicit behaviour and it might be better to
        //      just abort, as is currently the case.
        //      In order to implement the Jacobian assembly entirely, in addition
        //      to the operator() in VectorMatrixAssembler there has to be a method
        //      that dispatches the Jacobian assembly.
        //      Similarly, then the NeumannBC specialization of VectorMatrixAssembler
        //      probably can be merged into the main class s.t. one has only one
        //      type of VectorMatrixAssembler (for each equation type) with the
        //      three methods assemble(), assembleJacobian() and assembleNeumannBC().
        //      That list can be extended, e.g. by methods for the assembly of
        //      source terms.
        //      UPDATE: Probably it is better to keep a separate NeumannBC version of the
        //      VectoMatrixAssembler since that will work for all kinds of processes.
    }

    std::vector<DirichletBc<Index>> const* getKnownSolutions(
        double const /*t*/) const override final
    {
        return &_dirichlet_bcs;
    }

    NonlinearSolver& getNonlinearSolver() const
    {
        return _nonlinear_solver;
    }

    TimeDiscretization& getTimeDiscretization() const
    {
        return *_time_discretization;
    }

private:
    /// Process specific initialization called by initialize().
    virtual void initializeConcreteProcess(
        NumLib::LocalToGlobalIndexMap const& dof_table,
        MeshLib::Mesh const& mesh,
        unsigned const integration_order) = 0;

    virtual void assembleConcreteProcess(
        const double t, GlobalVector const& x,
        GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b) = 0;

    virtual void assembleJacobianConcreteProcess(
        const double /*t*/, GlobalVector const& /*x*/, GlobalVector const& /*xdot*/,
        const double /*dxdot_dx*/, GlobalMatrix const& /*M*/,
        const double /*dx_dx*/, GlobalMatrix const& /*K*/,
        GlobalMatrix& /*Jac*/)
    {
        OGS_FATAL("The concrete implementation of this Process did not override the"
            " assembleJacobianConcreteProcess() method."
            " Hence, no analytical Jacobian is provided for this process"
            " and the Newton-Raphson method cannot be used to solve it.");
    }

    void constructDofTable()
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

    /// Sets the initial condition values in the solution vector x for a given
    /// process variable and component.
    void setInitialConditions(ProcessVariable const& variable,
                              int const variable_id,
                              int const component_id,
                              GlobalVector& x)
    {
        std::size_t const n_nodes = _mesh.getNumberOfNodes();
        for (std::size_t node_id = 0; node_id < n_nodes; ++node_id)
        {
            MeshLib::Location const l(_mesh.getID(),
                                      MeshLib::MeshItemType::Node, node_id);
            auto global_index =
                std::abs(_local_to_global_index_map->getGlobalIndex(
                    l, variable_id, component_id));
#ifdef USE_PETSC
            // The global indices of the ghost entries of the global matrix or
            // the global vectors need to be set as negative values for equation
            // assembly, however the global indices start from zero.  Therefore,
            // any ghost entry with zero index is assigned an negative value of
            // the vector size or the matrix dimension.  To assign the initial
            // value for the ghost entries, the negative indices of the ghost
            // entries are restored to zero.
            if (global_index == x.size())
                global_index = 0;
#endif
            x.set(global_index,
                  variable.getInitialConditionValue(node_id, component_id));
        }
    }

    void createDirichletBcs(ProcessVariable& variable, int const variable_id,
                            int const component_id)
    {
        MeshGeoToolsLib::MeshNodeSearcher& mesh_node_searcher =
            MeshGeoToolsLib::MeshNodeSearcher::getMeshNodeSearcher(
                variable.getMesh());

        variable.initializeDirichletBCs(std::back_inserter(_dirichlet_bcs),
                                        mesh_node_searcher,
                                        *_local_to_global_index_map,
                                        variable_id,
                                        component_id);
    }

    void createNeumannBcs(ProcessVariable& variable, int const variable_id,
                          int const component_id)
    {
        // Find mesh nodes.
        MeshGeoToolsLib::MeshNodeSearcher& mesh_node_searcher =
            MeshGeoToolsLib::MeshNodeSearcher::getMeshNodeSearcher(
                variable.getMesh());
        MeshGeoToolsLib::BoundaryElementsSearcher mesh_element_searcher(
            variable.getMesh(), mesh_node_searcher);

        // Create a neumann BC for the process variable storing them in the
        // _neumann_bcs vector.
        variable.createNeumannBcs<GlobalSetup>(
                                  std::back_inserter(_neumann_bcs),
                                  mesh_element_searcher,
                                  _integration_order,
                                  *_local_to_global_index_map,
                                  variable_id,
                                  component_id);
    }

    /// Computes and stores global matrix' sparsity pattern from given
    /// DOF-table.
    void computeSparsityPattern()
    {
        _sparsity_pattern = std::move(NumLib::computeSparsityPattern(
            *_local_to_global_index_map, _mesh));
    }

protected:
    MeshLib::Mesh& _mesh;
    std::unique_ptr<MeshLib::MeshSubset const> _mesh_subset_all_nodes;

    std::unique_ptr<NumLib::LocalToGlobalIndexMap>
        _local_to_global_index_map;

    SecondaryVariableCollection<GlobalVector> _secondary_variables;
    ProcessOutput<GlobalVector> _process_output;

private:
    unsigned const _integration_order = 2;
    GlobalSparsityPattern _sparsity_pattern;

    std::vector<DirichletBc<GlobalIndexType>> _dirichlet_bcs;
    std::vector<std::unique_ptr<NeumannBc<GlobalSetup>>> _neumann_bcs;

    NonlinearSolver& _nonlinear_solver;
    std::unique_ptr<TimeDiscretization> _time_discretization;

    /// Variables used by this process.
    std::vector<std::reference_wrapper<ProcessVariable>> _process_variables;
};

/// Find process variables in \c variables whose names match the settings under
/// the given \c tag_names in the \c process_config.
///
/// In the process config a process variable is referenced by a name. For
/// example it will be looking for a variable named "H" in the list of process
/// variables when the tag is "hydraulic_head":
/// \code
///     <process>
///         ...
///         <process_variables>
///             <hydraulic_head>H</hydraulic_head>
///             ...
///         </process_variables>
///         ...
///     </process>
/// \endcode
///
/// \return a vector of references to the found variable(s).
std::vector<std::reference_wrapper<ProcessVariable>>
findProcessVariables(
        std::vector<ProcessVariable> const& variables,
        BaseLib::ConfigTree const& process_config,
        std::initializer_list<std::string> tag_names);

/// Find a parameter of specific type for a name given in the process
/// configuration under the tag.
/// In the process config a parameter is referenced by a name. For example it
/// will be looking for a parameter named "K" in the list of parameters
/// when the tag is "hydraulic_conductivity":
/// \code
///     <process>
///         ...
///         <hydraulic_conductivity>K</hydraulic_conductivity>
///     </process>
/// \endcode
/// and return a reference to that parameter. Additionally it checks for the
/// type of the found parameter.
template <typename... ParameterArgs>
Parameter<ParameterArgs...>& findParameter(
    BaseLib::ConfigTree const& process_config, std::string const& tag,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters)
{
    // Find parameter name in process config.
    //! \ogs_file_special
    auto const name = process_config.getConfigParameter<std::string>(tag);

    // Find corresponding parameter by name.
    auto const parameter_it =
        std::find_if(parameters.cbegin(), parameters.cend(),
                     [&name](std::unique_ptr<ParameterBase> const& p)
                     {
                         return p->name == name;
                     });

    if (parameter_it == parameters.end())
    {
        OGS_FATAL(
            "Could not find parameter '%s' in the provided parameters list for "
            "config tag <%s>.",
            name.c_str(), tag.c_str());
    }
    DBUG("Found parameter \'%s\'.", (*parameter_it)->name.c_str());

    // Check the type correctness of the found parameter.
    auto* const parameter =
        dynamic_cast<Parameter<ParameterArgs...>*>(parameter_it->get());
    if (!parameter)
    {
        OGS_FATAL("The read parameter is of incompatible type.");
    }
    return *parameter;
}

}  // namespace ProcessLib

#endif  // PROCESS_LIB_PROCESS_H_
