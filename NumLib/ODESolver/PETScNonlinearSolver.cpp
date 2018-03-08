/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#ifdef USE_PETSC

#include "PETScNonlinearSolver.h"

#include <petscmat.h>
#include <petscvec.h>

#include "BaseLib/RunTime.h"

namespace
{
struct PetscContext
{
    using System = NumLib::NonlinearSystem<NumLib::NonlinearSolverTag::Newton>;
    System* system;
    std::vector<GlobalVector*>& x;
    std::vector<GlobalVector*> const& x_prev;
    GlobalVector* r;
    GlobalMatrix* J;
    int const process_id;
};

PetscErrorCode updateResidual(SNES /*snes*/, Vec x, Vec petsc_r,
                              void* petsc_context)
{
    auto context = static_cast<PetscContext*>(petsc_context);

    DBUG("PETScNonlinearSolver: residual callback called.");

    VecCopy(x, context->x[context->process_id]->getRawVector());

    // Assemble in ogs context.
    BaseLib::RunTime time_assembly;
    time_assembly.start();
    context->system->assemble(context->x, context->x_prev, context->process_id);

    INFO("[time] Assembly took {} s.", time_assembly.elapsed());
    context->system->getResidual(*context->x[context->process_id],
                                 *context->x_prev[context->process_id],
                                 *context->r);
    context->J->finalizeAssembly();

    context->system->getJacobian(*context->J);
    context->system->applyKnownSolutionsNewton(
        *context->J, *context->r, *context->x[context->process_id]);

    VecCopy(context->r->getRawVector(), petsc_r);

    return 0;
}

PetscErrorCode updateJacobian(SNES /*snes*/, Vec /*x*/, Mat J,
                              Mat /*same as J*/, void* petsc_context)
{
    DBUG("PETScNonlinearSolver: jacobian callback called.");
    // Assume the system is already assembled.
    // Copy the results into petsc vectors.

    auto context = static_cast<PetscContext*>(petsc_context);
    MatAssemblyBegin(J, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(J, MAT_FINAL_ASSEMBLY);
    MatCopy(context->J->getRawMatrix(), J, DIFFERENT_NONZERO_PATTERN);
    return 0;
}
}  // namespace

namespace NumLib
{
PETScNonlinearSolver::PETScNonlinearSolver(
    GlobalLinearSolver& /*linear_solver*/)
{
    SNESCreate(PETSC_COMM_WORLD, &_snes_solver);
    SNESSetFromOptions(_snes_solver);
#ifndef NDEBUG
    PetscOptionsView(nullptr, PETSC_VIEWER_STDOUT_WORLD);
#endif  // NDEBUG
}

void PETScNonlinearSolver::setEquationSystem(System& eq,
                                             ConvergenceCriterion& conv_crit)
{
    _equation_system = &eq;
    _convergence_criterion = &conv_crit;
}

void PETScNonlinearSolver::compensateNonEquilibriumInitialResiduum(
    bool const value)
{
    _compensate_non_equilibrium_initial_residuum = value;
}

void PETScNonlinearSolver::calculateNonEquilibriumInitialResiduum(
    std::vector<GlobalVector*> const& /*x*/,
    std::vector<GlobalVector*> const& /*x_prev*/, int const /*process_id*/)
{
    if (!_compensate_non_equilibrium_initial_residuum)
    {
        return;
    }

    OGS_FATAL(
        "Non-equilibrium initial residuum is not implemented for the "
        "PETScNonlinearSolver.");
}

NonlinearSolverStatus PETScNonlinearSolver::solve(
    std::vector<GlobalVector*>& x,
    std::vector<GlobalVector*> const& x_prev,
    std::function<void(
        int,
        std::vector<GlobalVector*> const&)> const& /*postIterationCallback*/,
    int const process_id)
{
    DBUG("PETScNonlinearSolver: solve()");
    using TimeDiscretizedSystem =
        TimeDiscretizedODESystem<ODESystemTag::FirstOrderImplicitQuasilinear,
                                 NonlinearSolverTag::Newton>;

    auto* system = static_cast<TimeDiscretizedSystem*>(_equation_system);

    DBUG("PETScNonlinearSolver: create vectors");
    // r and J on which the ogs assembly operates.
    auto& J = NumLib::GlobalMatrixProvider::provider.getMatrix(
        system->getMatrixSpecifications(process_id), _jacobian_id);
    auto& r = NumLib::GlobalVectorProvider::provider.getVector(
        system->getMatrixSpecifications(process_id), _residual_id);

    // Vectors and matrices used by SNES for solutions. These will be locked
    // during the SNES' solve call.
    auto& J_snes = NumLib::GlobalMatrixProvider::provider.getMatrix(
        system->getMatrixSpecifications(process_id), _petsc_jacobian_id);
    auto& r_snes = NumLib::GlobalVectorProvider::provider.getVector(
        system->getMatrixSpecifications(process_id), _petsc_residual_id);
    auto& x_snes = NumLib::GlobalVectorProvider::provider.getVector(
        system->getMatrixSpecifications(process_id), _petsc_x_id);
    MathLib::LinAlg::copy(*x[process_id], x_snes);  // Initial guess.

    PetscContext petsc_context{_equation_system, x, x_prev, &r, &J, process_id};

    DBUG("PETScNonlinearSolver: set function");
    SNESSetFunction(_snes_solver, r_snes.getRawVector(), updateResidual,
                    &petsc_context);

    DBUG("PETScNonlinearSolver: set jacobian");
    // The jacobian and the preconditioner matrices point to the same location.
    SNESSetJacobian(_snes_solver, J_snes.getRawMatrix(), J_snes.getRawMatrix(),
                    updateJacobian, &petsc_context);

    std::unique_ptr<GlobalVector> xl = nullptr;
    std::unique_ptr<GlobalVector> xu = nullptr;

    SNESType snes_type;
    SNESGetType(_snes_solver, &snes_type);
    if ((std::strcmp(snes_type, SNESVINEWTONRSLS) == 0) ||
        (std::strcmp(snes_type, SNESVINEWTONSSLS) == 0))
    {
        // Set optional constraints via callback.
        DBUG("PETScNonlinearSolver: set constraints");
        xl = MathLib::MatrixVectorTraits<GlobalVector>::newInstance(
            system->getMatrixSpecifications(process_id));
        xu = MathLib::MatrixVectorTraits<GlobalVector>::newInstance(
            system->getMatrixSpecifications(process_id));

        system->updateConstraints(*xl, *xu, process_id);
        MathLib::finalizeVectorAssembly(*xl);
        MathLib::finalizeVectorAssembly(*xu);

        SNESVISetVariableBounds(_snes_solver, xl->getRawVector(),
                                xu->getRawVector());
    }

    DBUG("PETScNonlinearSolver: call SNESSolve");
    SNESSolve(_snes_solver, nullptr, x_snes.getRawVector());

    SNESConvergedReason reason = SNES_CONVERGED_ITERATING;
    SNESGetConvergedReason(_snes_solver, &reason);
    INFO("PETSsSNES convergence reason {}.", reason);

    PetscInt iterations;
    SNESGetIterationNumber(_snes_solver, &iterations);
    INFO("PETScSNES used {} iterations.", iterations);

    // Copy back the solution.
    MathLib::LinAlg::copy(x_snes, *x[process_id]);

    NumLib::GlobalVectorProvider::provider.releaseVector(x_snes);
    NumLib::GlobalVectorProvider::provider.releaseVector(r_snes);
    NumLib::GlobalMatrixProvider::provider.releaseMatrix(J_snes);
    NumLib::GlobalVectorProvider::provider.releaseVector(r);
    NumLib::GlobalMatrixProvider::provider.releaseMatrix(J);

    return {reason >= 0, iterations};
}

}  // namespace NumLib
#endif  // USE_PETSC
