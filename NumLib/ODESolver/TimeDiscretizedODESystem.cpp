/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "TimeDiscretizedODESystem.h"

#include "MathLib/LinAlg/ApplyKnownSolution.h"
#include "MathLib/LinAlg/UnifiedMatrixSetters.h"
#include "NumLib/IndexValueVector.h"
#include "NumLib/Exceptions.h"

namespace detail
{
//! Applies known solutions to the solution vector \c x.
template <typename Solutions, typename Vector>
void applyKnownSolutions(std::vector<Solutions> const* const known_solutions,
                         Vector& x)
{
    if (!known_solutions)
    {
        return;
    }

    for (auto const& bc : *known_solutions)
    {
        for (std::size_t i = 0; i < bc.ids.size(); ++i)
        {
            // TODO that might have bad performance for some Vector types, e.g.,
            // PETSc.
            MathLib::setVector(x, bc.ids[i], bc.values[i]);
        }
    }
    MathLib::LinAlg::finalizeAssembly(x);
}
}  // namespace detail

namespace NumLib
{
TimeDiscretizedODESystem<ODESystemTag::FirstOrderImplicitQuasilinear,
                         NonlinearSolverTag::Newton>::
    TimeDiscretizedODESystem(const int process_id, ODE& ode,
                             TimeDisc& time_discretization)
    : ode_(ode),
      time_disc_(time_discretization),
      mat_trans_(createMatrixTranslator<ODETag>(time_discretization))
{
    Jac_ = &NumLib::GlobalMatrixProvider::provider.getMatrix(
        ode_.getMatrixSpecifications(process_id), Jac_id_);
    M_ = &NumLib::GlobalMatrixProvider::provider.getMatrix(
        ode_.getMatrixSpecifications(process_id), M_id_);
    K_ = &NumLib::GlobalMatrixProvider::provider.getMatrix(
        ode_.getMatrixSpecifications(process_id), K_id_);
    b_ = &NumLib::GlobalVectorProvider::provider.getVector(
        ode_.getMatrixSpecifications(process_id), b_id_);
}

TimeDiscretizedODESystem<
    ODESystemTag::FirstOrderImplicitQuasilinear,
    NonlinearSolverTag::Newton>::~TimeDiscretizedODESystem()
{
    NumLib::GlobalMatrixProvider::provider.releaseMatrix(*Jac_);
    NumLib::GlobalMatrixProvider::provider.releaseMatrix(*M_);
    NumLib::GlobalMatrixProvider::provider.releaseMatrix(*K_);
    NumLib::GlobalVectorProvider::provider.releaseVector(*b_);
}

void TimeDiscretizedODESystem<ODESystemTag::FirstOrderImplicitQuasilinear,
                              NonlinearSolverTag::Newton>::
    assemble(std::vector<GlobalVector*> const& x_new_timestep,
             std::vector<GlobalVector*> const& x_prev,
             int const process_id)
{
    namespace LinAlg = MathLib::LinAlg;

    auto const t = time_disc_.getCurrentTime();
    auto const dt = time_disc_.getCurrentTimeIncrement();
    auto const& x_curr = *x_new_timestep[process_id];
    auto const dxdot_dx = time_disc_.getNewXWeight();

    std::vector<GlobalVector*> xdot(x_new_timestep.size());
    for (auto& v : xdot)
    {
        v = &NumLib::GlobalVectorProvider::provider.getVector();
        time_disc_.getXdot(*x_new_timestep[process_id], *x_prev[process_id],
                           *v);
    }

    M_->setZero();
    K_->setZero();
    b_->setZero();
    Jac_->setZero();

    ode_.preAssemble(t, dt, x_curr);
    try
    {
        ode_.assembleWithJacobian(t, dt, x_new_timestep, *xdot[process_id],
                                  dxdot_dx, 1.0, process_id, *M_, *K_, *b_,
                                  *Jac_);
    }
    catch (AssemblyException const&)
    {
        for (auto& v : xdot)
        {
            NumLib::GlobalVectorProvider::provider.releaseVector(*v);
        }
        throw;
    }

    LinAlg::finalizeAssembly(*M_);
    LinAlg::finalizeAssembly(*K_);
    LinAlg::finalizeAssembly(*b_);
    MathLib::LinAlg::finalizeAssembly(*Jac_);

    for (auto& v : xdot)
    {
        NumLib::GlobalVectorProvider::provider.releaseVector(*v);
    }
}

void TimeDiscretizedODESystem<
    ODESystemTag::FirstOrderImplicitQuasilinear,
    NonlinearSolverTag::Newton>::getResidual(GlobalVector const& x_new_timestep,
                                             GlobalVector const& x_prev,
                                             GlobalVector& res) const
{
    // TODO Maybe the duplicate calculation of xdot here and in assembleJacobian
    //      can be optimuized. However, that would make the interface a bit more
    //      fragile.
    auto& xdot = NumLib::GlobalVectorProvider::provider.getVector(xdot_id_);
    time_disc_.getXdot(x_new_timestep, x_prev, xdot);

    mat_trans_->computeResidual(*M_, *K_, *b_, x_new_timestep, xdot, res);

    NumLib::GlobalVectorProvider::provider.releaseVector(xdot);
}

void TimeDiscretizedODESystem<
    ODESystemTag::FirstOrderImplicitQuasilinear,
    NonlinearSolverTag::Newton>::getJacobian(GlobalMatrix& Jac) const
{
    mat_trans_->computeJacobian(*Jac_, Jac);
}

void TimeDiscretizedODESystem<
    ODESystemTag::FirstOrderImplicitQuasilinear,
    NonlinearSolverTag::Newton>::computeKnownSolutions(GlobalVector const& x,
                                                       int const process_id)
{
    known_solutions_ =
        ode_.getKnownSolutions(time_disc_.getCurrentTime(), x, process_id);
}

void TimeDiscretizedODESystem<
    ODESystemTag::FirstOrderImplicitQuasilinear,
    NonlinearSolverTag::Newton>::applyKnownSolutions(GlobalVector& x) const
{
    ::detail::applyKnownSolutions(known_solutions_, x);
}

void TimeDiscretizedODESystem<ODESystemTag::FirstOrderImplicitQuasilinear,
                              NonlinearSolverTag::Newton>::
    applyKnownSolutionsNewton(GlobalMatrix& Jac, GlobalVector& res,
                              GlobalVector& minus_delta_x) const
{
    if (!known_solutions_)
    {
        return;
    }

    using IndexType = MathLib::MatrixVectorTraits<GlobalMatrix>::Index;
    std::vector<IndexType> ids;
    for (auto const& bc : *known_solutions_)
    {
        std::copy(bc.ids.cbegin(), bc.ids.cend(), std::back_inserter(ids));
    }

    // For the Newton method the values must be zero
    std::vector<double> values(ids.size(), 0);
    MathLib::applyKnownSolution(Jac, res, minus_delta_x, ids, values);
}

TimeDiscretizedODESystem<ODESystemTag::FirstOrderImplicitQuasilinear,
                         NonlinearSolverTag::Picard>::
    TimeDiscretizedODESystem(const int process_id, ODE& ode,
                             TimeDisc& time_discretization)
    : ode_(ode),
      time_disc_(time_discretization),
      mat_trans_(createMatrixTranslator<ODETag>(time_discretization))
{
    M_ = &NumLib::GlobalMatrixProvider::provider.getMatrix(
        ode.getMatrixSpecifications(process_id), M_id_);
    K_ = &NumLib::GlobalMatrixProvider::provider.getMatrix(
        ode.getMatrixSpecifications(process_id), K_id_);
    b_ = &NumLib::GlobalVectorProvider::provider.getVector(
        ode.getMatrixSpecifications(process_id), b_id_);
}

TimeDiscretizedODESystem<
    ODESystemTag::FirstOrderImplicitQuasilinear,
    NonlinearSolverTag::Picard>::~TimeDiscretizedODESystem()
{
    NumLib::GlobalMatrixProvider::provider.releaseMatrix(*M_);
    NumLib::GlobalMatrixProvider::provider.releaseMatrix(*K_);
    NumLib::GlobalVectorProvider::provider.releaseVector(*b_);
}

void TimeDiscretizedODESystem<ODESystemTag::FirstOrderImplicitQuasilinear,
                              NonlinearSolverTag::Picard>::
    assemble(std::vector<GlobalVector*> const& x_new_timestep,
             std::vector<GlobalVector*> const& x_prev,
             int const process_id)
{
    namespace LinAlg = MathLib::LinAlg;

    auto const t = time_disc_.getCurrentTime();
    auto const dt = time_disc_.getCurrentTimeIncrement();
    auto const& x_curr = *x_new_timestep[process_id];
    std::vector<GlobalVector*> xdot(x_new_timestep.size());
    for (auto& v : xdot)
    {
        v = &NumLib::GlobalVectorProvider::provider.getVector();
        time_disc_.getXdot(*x_new_timestep[process_id], *x_prev[process_id],
                           *v);
    }

    M_->setZero();
    K_->setZero();
    b_->setZero();

    ode_.preAssemble(t, dt, x_curr);
    ode_.assemble(t, dt, x_new_timestep, xdot, process_id, *M_, *K_, *b_);

    LinAlg::finalizeAssembly(*M_);
    LinAlg::finalizeAssembly(*K_);
    LinAlg::finalizeAssembly(*b_);
}

void TimeDiscretizedODESystem<
    ODESystemTag::FirstOrderImplicitQuasilinear,
    NonlinearSolverTag::Picard>::computeKnownSolutions(GlobalVector const& x,
                                                       int const process_id)
{
    known_solutions_ =
        ode_.getKnownSolutions(time_disc_.getCurrentTime(), x, process_id);
}

void TimeDiscretizedODESystem<
    ODESystemTag::FirstOrderImplicitQuasilinear,
    NonlinearSolverTag::Picard>::applyKnownSolutions(GlobalVector& x) const
{
    ::detail::applyKnownSolutions(known_solutions_, x);
}

void TimeDiscretizedODESystem<ODESystemTag::FirstOrderImplicitQuasilinear,
                              NonlinearSolverTag::Picard>::
    applyKnownSolutionsPicard(GlobalMatrix& A,
                              GlobalVector& rhs,
                              GlobalVector& x) const
{
    if (!known_solutions_)
    {
        return;
    }

    using IndexType = MathLib::MatrixVectorTraits<GlobalMatrix>::Index;
    std::vector<IndexType> ids;
    std::vector<double> values;
    for (auto const& bc : *known_solutions_)
    {
        std::copy(bc.ids.cbegin(), bc.ids.cend(), std::back_inserter(ids));
        std::copy(bc.values.cbegin(), bc.values.cend(),
                  std::back_inserter(values));
    }
    MathLib::applyKnownSolution(A, rhs, x, ids, values);
}

}  // namespace NumLib
