/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "TimeDiscretizedODESystem.h"

#include <range/v3/numeric/accumulate.hpp>
#include <range/v3/view/transform.hpp>

#include "MathLib/LinAlg/ApplyKnownSolution.h"
#include "MathLib/LinAlg/UnifiedMatrixSetters.h"
#include "NumLib/Exceptions.h"
#include "NumLib/IndexValueVector.h"

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
    _b = &NumLib::GlobalVectorProvider::provider.getVector(
        _ode.getMatrixSpecifications(process_id), _b_id);
}

TimeDiscretizedODESystem<
    ODESystemTag::FirstOrderImplicitQuasilinear,
    NonlinearSolverTag::Newton>::~TimeDiscretizedODESystem()
{
    NumLib::GlobalMatrixProvider::provider.releaseMatrix(*_Jac);
    NumLib::GlobalVectorProvider::provider.releaseVector(*_b);
}

void TimeDiscretizedODESystem<ODESystemTag::FirstOrderImplicitQuasilinear,
                              NonlinearSolverTag::Newton>::
    assemble(std::vector<GlobalVector*> const& x_new_timestep,
             std::vector<GlobalVector*> const& x_prev,
             int const process_id)
{
    auto const t = _time_disc.getCurrentTime();
    auto const dt = _time_disc.getCurrentTimeIncrement();
    auto const& x_curr = *x_new_timestep[process_id];

    _b->setZero();
    _Jac->setZero();

    _ode.preAssemble(t, dt, x_curr);
    _ode.assembleWithJacobian(t, dt, x_new_timestep, x_prev, process_id, *_b,
                              *_Jac);

    MathLib::LinAlg::finalizeAssembly(*_b);
    MathLib::LinAlg::finalizeAssembly(*_Jac);
}

void TimeDiscretizedODESystem<ODESystemTag::FirstOrderImplicitQuasilinear,
                              NonlinearSolverTag::Newton>::
    getResidual(GlobalVector const& /*x_new_timestep*/,
                GlobalVector const& /*x_prev*/,
                GlobalVector& res) const
{
    MathLib::LinAlg::copy(*_b, res);   // res = b
    MathLib::LinAlg::scale(res, -1.);  // res = -b
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
                              GlobalVector const& x,
                              GlobalVector& minus_delta_x) const
{
    if (!_known_solutions)
    {
        return;
    }

    using IndexType = MathLib::MatrixVectorTraits<GlobalMatrix>::Index;
    std::size_t const size = ranges::accumulate(
        *_known_solutions | ranges::views::transform([](auto const& bc)
                                                     { return bc.ids.size(); }),
        0);
    std::vector<IndexType> ids;
    ids.reserve(size);
    std::vector<double> values;
    values.reserve(size);

    for (auto const& bc : *_known_solutions)
    {
        for (std::size_t i = 0; i < bc.ids.size(); ++i)
        {
            auto const id = bc.ids[i];
            ids.push_back(id);
            // minus_delta_x will be set to the difference between the current
            // value and the Dirichlet BC value.
            values.push_back(x[id] - bc.values[i]);
        }
    }

    MathLib::applyKnownSolution(Jac, res, minus_delta_x, ids, values);
}

void TimeDiscretizedODESystem<ODESystemTag::FirstOrderImplicitQuasilinear,
                              NonlinearSolverTag::Newton>::
    applyKnownSolutionsPETScSNES(GlobalMatrix& Jac, GlobalVector& res,
                                 GlobalVector& x) const
{
    if (!_known_solutions)
    {
        return;
    }

    using IndexType = MathLib::MatrixVectorTraits<GlobalMatrix>::Index;
    std::vector<IndexType> ids;
    for (auto const& bc : *_known_solutions)
    {
        ids.insert(end(ids), begin(bc.ids), end(bc.ids));
    }

    // For the Newton method the values must be zero
    std::vector<double> values(ids.size(), 0);
    MathLib::applyKnownSolution(Jac, res, x, ids, values);
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

    _M->setZero();
    _K->setZero();
    _b->setZero();

    _ode.preAssemble(t, dt, x_curr);
    _ode.assemble(t, dt, x_new_timestep, x_prev, process_id, *_M, *_K, *_b);

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
        ids.insert(end(ids), begin(bc.ids), end(bc.ids));
        values.insert(end(values), begin(bc.values), end(bc.values));
    }
    MathLib::applyKnownSolution(A, rhs, x, ids, values);
}

}  // namespace NumLib
