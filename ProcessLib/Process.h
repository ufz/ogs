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

#include "NumLib/ODESolver/NonlinearSolver.h"
#include "NumLib/ODESolver/ODESystem.h"
#include "NumLib/ODESolver/TimeDiscretization.h"
#include "NumLib/NamedFunctionCaller.h"
#include "ProcessLib/BoundaryCondition/BoundaryConditionCollection.h"

#include "ExtrapolatorData.h"
#include "Parameter.h"
#include "ProcessOutput.h"
#include "SecondaryVariable.h"
#include "CachedSecondaryVariable.h"

namespace MeshLib
{
class Mesh;
}

namespace ProcessLib
{
class Process
    : public NumLib::ODESystem<  // TODO: later on use a simpler ODE system
          NumLib::ODESystemTag::FirstOrderImplicitQuasilinear,
          NumLib::NonlinearSolverTag::Newton>
{
public:
    using NonlinearSolver = NumLib::NonlinearSolverBase;
    using TimeDiscretization = NumLib::TimeDiscretization;

    Process(
        MeshLib::Mesh& mesh,
        NonlinearSolver& nonlinear_solver,
        std::unique_ptr<TimeDiscretization>&& time_discretization,
        std::unique_ptr<NumLib::ConvergenceCriterion>&& convergence_criterion,
        std::vector<std::reference_wrapper<ProcessVariable>>&&
            process_variables,
        SecondaryVariableCollection&& secondary_variables,
        ProcessOutput&& process_output,
        NumLib::NamedFunctionCaller&& named_function_caller);

    /// Preprocessing before starting assembly for new timestep.
    virtual void preTimestep(GlobalVector const& /*x*/, const double /*t*/,
                             const double /*delta_t*/)
    {
    }

    /// Postprocessing after a complete timestep.
    virtual void postTimestep(GlobalVector const& /*x*/) {}

    void preIteration(const unsigned iter,
                      GlobalVector const& x) override final;

    NumLib::IterationResult postIteration(GlobalVector const& x) override final;

    /// Process output.
    /// The file_name is indicating the name of possible output file.
    void output(std::string const& file_name,
                const unsigned /*timestep*/,
                GlobalVector const& x) const;

    void initialize();

    void setInitialConditions(const double t, GlobalVector& x);

    MathLib::MatrixSpecifications getMatrixSpecifications()
        const override final;

    void assemble(const double t, GlobalVector const& x, GlobalMatrix& M,
                  GlobalMatrix& K, GlobalVector& b) override final;

    void assembleJacobian(const double t, GlobalVector const& x,
                          GlobalVector const& xdot, const double dxdot_dx,
                          GlobalMatrix const& M, const double dx_dx,
                          GlobalMatrix const& K,
                          GlobalMatrix& Jac) override final;
    std::vector<NumLib::IndexValueVector<GlobalIndexType>> const*
    getKnownSolutions(

        double const t) const override final
    {
        return _boundary_conditions.getKnownSolutions(t);
    }

    NonlinearSolver& getNonlinearSolver() const { return _nonlinear_solver; }
    TimeDiscretization& getTimeDiscretization() const
    {
        return *_time_discretization;
    }

    NumLib::ConvergenceCriterion& getConvergenceCriterion() const
    {
        return *_convergence_criterion;
    }

protected:
    NumLib::Extrapolator& getExtrapolator() const
    {
        return _extrapolator_data.getExtrapolator();
    }

    NumLib::LocalToGlobalIndexMap const& getSingleComponentDOFTable() const
    {
        return _extrapolator_data.getDOFTable();
    }

private:
    /// Process specific initialization called by initialize().
    virtual void initializeConcreteProcess(
        NumLib::LocalToGlobalIndexMap const& dof_table,
        MeshLib::Mesh const& mesh,
        unsigned const integration_order) = 0;

    virtual void assembleConcreteProcess(const double t, GlobalVector const& x,
                                         GlobalMatrix& M, GlobalMatrix& K,
                                         GlobalVector& b) = 0;

    virtual void assembleJacobianConcreteProcess(
        const double /*t*/, GlobalVector const& /*x*/,
        GlobalVector const& /*xdot*/, const double /*dxdot_dx*/,
        GlobalMatrix const& /*M*/, const double /*dx_dx*/,
        GlobalMatrix const& /*K*/, GlobalMatrix& /*Jac*/);

    virtual void preIterationConcreteProcess(const unsigned /*iter*/,
                                             GlobalVector const& /*x*/){}

    virtual NumLib::IterationResult postIterationConcreteProcess(
        GlobalVector const& /*x*/)
    {
        return NumLib::IterationResult::SUCCESS;
    }

    void constructDofTable();

    void initializeExtrapolator();

    /// Finishes the \c _named_function_caller and adds a secondary variable for
    /// each of the named functions.
    void finishNamedFunctionsInitialization();

    /// Computes and stores global matrix' sparsity pattern from given
    /// DOF-table.
    void computeSparsityPattern();

protected:
    MeshLib::Mesh& _mesh;
    std::unique_ptr<MeshLib::MeshSubset const> _mesh_subset_all_nodes;

    std::unique_ptr<NumLib::LocalToGlobalIndexMap> _local_to_global_index_map;

    SecondaryVariableCollection _secondary_variables;
    ProcessOutput _process_output;

    NumLib::NamedFunctionCaller _named_function_caller;
    std::vector<std::unique_ptr<CachedSecondaryVariable>>
        _cached_secondary_variables;
    SecondaryVariableContext _secondary_variable_context;

private:
    unsigned const _integration_order = 2;
    GlobalSparsityPattern _sparsity_pattern;

    NonlinearSolver& _nonlinear_solver;
    std::unique_ptr<TimeDiscretization> _time_discretization;
    std::unique_ptr<NumLib::ConvergenceCriterion> _convergence_criterion;

    /// Variables used by this process.
    std::vector<std::reference_wrapper<ProcessVariable>> _process_variables;

    BoundaryConditionCollection _boundary_conditions;

    ExtrapolatorData _extrapolator_data;
};

}  // namespace ProcessLib

#endif  // PROCESS_LIB_PROCESS_H_
