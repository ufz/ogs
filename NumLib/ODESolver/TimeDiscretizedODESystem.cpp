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
    : _ode(ode),
      _time_disc(time_discretization),
      _mat_trans(createMatrixTranslator<ODETag>(time_discretization))
{
    _Jac = &NumLib::GlobalMatrixProvider::provider.getMatrix(
        _ode.getMatrixSpecifications(process_id), _Jac_id);
    _M = &NumLib::GlobalMatrixProvider::provider.getMatrix(
        _ode.getMatrixSpecifications(process_id), _M_id);
    _K = &NumLib::GlobalMatrixProvider::provider.getMatrix(
        _ode.getMatrixSpecifications(process_id), _K_id);
    _b = &NumLib::GlobalVectorProvider::provider.getVector(
        _ode.getMatrixSpecifications(process_id), _b_id);
}

TimeDiscretizedODESystem<
    ODESystemTag::FirstOrderImplicitQuasilinear,
    NonlinearSolverTag::Newton>::~TimeDiscretizedODESystem()
{
    NumLib::GlobalMatrixProvider::provider.releaseMatrix(*_Jac);
    NumLib::GlobalMatrixProvider::provider.releaseMatrix(*_M);
    NumLib::GlobalMatrixProvider::provider.releaseMatrix(*_K);
    NumLib::GlobalVectorProvider::provider.releaseVector(*_b);
}

void TimeDiscretizedODESystem<ODESystemTag::FirstOrderImplicitQuasilinear,
                              NonlinearSolverTag::Newton>::
    assemble(std::vector<GlobalVector*> const& x_new_timestep,
             std::vector<GlobalVector*> const& x_prev,
             int const process_id)
{
    namespace LinAlg = MathLib::LinAlg;

    auto const t = _time_disc.getCurrentTime();
    auto const dt = _time_disc.getCurrentTimeIncrement();
    auto const& x_curr = *x_new_timestep[process_id];
    auto const dxdot_dx = _time_disc.getNewXWeight();

    std::vector<GlobalVector*> xdot(x_new_timestep.size());
    for (std::size_t i = 0; i < xdot.size(); i++)
    {
        xdot[i] = &NumLib::GlobalVectorProvider::provider.getVector();
        _time_disc.getXdot(*x_new_timestep[i], *x_prev[i], *xdot[i]);
    }

    _M->setZero();
    _K->setZero();
    _b->setZero();
    _Jac->setZero();

    _ode.preAssemble(t, dt, x_curr);
    try
    {
        _ode.assembleWithJacobian(t, dt, x_new_timestep, xdot, dxdot_dx, 1.0,
                                  process_id, *_M, *_K, *_b, *_Jac);
    }
    catch (AssemblyException const&)
    {
        for (auto& v : xdot)
        {
            NumLib::GlobalVectorProvider::provider.releaseVector(*v);
        }
        throw;
    }

    LinAlg::finalizeAssembly(*_M);
    LinAlg::finalizeAssembly(*_K);
    LinAlg::finalizeAssembly(*_b);
    MathLib::LinAlg::finalizeAssembly(*_Jac);

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
    auto& xdot = NumLib::GlobalVectorProvider::provider.getVector(_xdot_id);
    _time_disc.getXdot(x_new_timestep, x_prev, xdot);

    _mat_trans->computeResidual(*_M, *_K, *_b, x_new_timestep, xdot, res);

    NumLib::GlobalVectorProvider::provider.releaseVector(xdot);
}

void TimeDiscretizedODESystem<
    ODESystemTag::FirstOrderImplicitQuasilinear,
    NonlinearSolverTag::Newton>::getJacobian(GlobalMatrix& Jac) const
{
    _mat_trans->computeJacobian(*_Jac, Jac);
}

void TimeDiscretizedODESystem<
    ODESystemTag::FirstOrderImplicitQuasilinear,
    NonlinearSolverTag::Newton>::computeKnownSolutions(GlobalVector const& x,
                                                       int const process_id)
{
    _known_solutions =
        _ode.getKnownSolutions(_time_disc.getCurrentTime(), x, process_id);
}

void TimeDiscretizedODESystem<
    ODESystemTag::FirstOrderImplicitQuasilinear,
    NonlinearSolverTag::Newton>::applyKnownSolutions(GlobalVector& x) const
{
    ::detail::applyKnownSolutions(_known_solutions, x);
}

void TimeDiscretizedODESystem<ODESystemTag::FirstOrderImplicitQuasilinear,
                              NonlinearSolverTag::Newton>::
    applyKnownSolutionsNewton(GlobalMatrix& Jac, GlobalVector& res,
                              GlobalVector& minus_delta_x) const
{
    if (!_known_solutions)
    {
        return;
    }

    using IndexType = MathLib::MatrixVectorTraits<GlobalMatrix>::Index;
    std::vector<IndexType> ids;
    for (auto const& bc : *_known_solutions)
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
    : _ode(ode),
      _time_disc(time_discretization),
      _mat_trans(createMatrixTranslator<ODETag>(time_discretization))
{
    _M = &NumLib::GlobalMatrixProvider::provider.getMatrix(
        ode.getMatrixSpecifications(process_id), _M_id);
    _K = &NumLib::GlobalMatrixProvider::provider.getMatrix(
        ode.getMatrixSpecifications(process_id), _K_id);
    _b = &NumLib::GlobalVectorProvider::provider.getVector(
        ode.getMatrixSpecifications(process_id), _b_id);
}

TimeDiscretizedODESystem<
    ODESystemTag::FirstOrderImplicitQuasilinear,
    NonlinearSolverTag::Picard>::~TimeDiscretizedODESystem()
{
    NumLib::GlobalMatrixProvider::provider.releaseMatrix(*_M);
    NumLib::GlobalMatrixProvider::provider.releaseMatrix(*_K);
    NumLib::GlobalVectorProvider::provider.releaseVector(*_b);
}

void TimeDiscretizedODESystem<ODESystemTag::FirstOrderImplicitQuasilinear,
                              NonlinearSolverTag::Picard>::
    assemble(std::vector<GlobalVector*> const& x_new_timestep,
             std::vector<GlobalVector*> const& x_prev,
             int const process_id)
{
    namespace LinAlg = MathLib::LinAlg;

    auto const t = _time_disc.getCurrentTime();
    auto const dt = _time_disc.getCurrentTimeIncrement();
    auto const& x_curr = *x_new_timestep[process_id];
    std::vector<GlobalVector*> xdot(x_new_timestep.size());
    for (std::size_t i = 0; i < xdot.size(); i++)
    {
        xdot[i] = &NumLib::GlobalVectorProvider::provider.getVector();
        _time_disc.getXdot(*x_new_timestep[i], *x_prev[i], *xdot[i]);
    }

    _M->setZero();
    _K->setZero();
    _b->setZero();

    _ode.preAssemble(t, dt, x_curr);
    _ode.assemble(t, dt, x_new_timestep, xdot, process_id, *_M, *_K, *_b);

    LinAlg::finalizeAssembly(*_M);
    LinAlg::finalizeAssembly(*_K);
    LinAlg::finalizeAssembly(*_b);
}

void TimeDiscretizedODESystem<
    ODESystemTag::FirstOrderImplicitQuasilinear,
    NonlinearSolverTag::Picard>::computeKnownSolutions(GlobalVector const& x,
                                                       int const process_id)
{
    _known_solutions =
        _ode.getKnownSolutions(_time_disc.getCurrentTime(), x, process_id);
}

void TimeDiscretizedODESystem<
    ODESystemTag::FirstOrderImplicitQuasilinear,
    NonlinearSolverTag::Picard>::applyKnownSolutions(GlobalVector& x) const
{
    ::detail::applyKnownSolutions(_known_solutions, x);
}

void TimeDiscretizedODESystem<ODESystemTag::FirstOrderImplicitQuasilinear,
                              NonlinearSolverTag::Picard>::
    applyKnownSolutionsPicard(GlobalMatrix& A,
                              GlobalVector& rhs,
                              GlobalVector& x) const
{
    if (!_known_solutions)
    {
        return;
    }

    using IndexType = MathLib::MatrixVectorTraits<GlobalMatrix>::Index;
    std::vector<IndexType> ids;
    std::vector<double> values;
    for (auto const& bc : *_known_solutions)
    {
        std::copy(bc.ids.cbegin(), bc.ids.cend(), std::back_inserter(ids));
        std::copy(bc.values.cbegin(), bc.values.cend(),
                  std::back_inserter(values));
    }
    MathLib::applyKnownSolution(A, rhs, x, ids, values);
}

}  // namespace NumLib
