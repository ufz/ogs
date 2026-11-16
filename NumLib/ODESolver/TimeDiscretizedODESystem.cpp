// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "TimeDiscretizedODESystem.h"

#include <range/v3/numeric/accumulate.hpp>
#include <range/v3/view/transform.hpp>

#include "BaseLib/Error.h"
#include "MathLib/LinAlg/ApplyKnownSolution.h"
#include "MathLib/LinAlg/UnifiedMatrixSetters.h"
#include "NumLib/DOF/GlobalMatrixProviders.h"
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
        MathLib::setVector(x, bc.ids, bc.values);
    }
    MathLib::LinAlg::finalizeAssembly(x);
}

//! Acquires a global matrix from the provider, which takes it back when the
//! last owner releases it.
std::shared_ptr<GlobalMatrix> provideMatrix(
    MathLib::MatrixSpecifications const& matrix_specifications)
{
    std::size_t id = 0u;
    return {&NumLib::GlobalMatrixProvider::provider.getMatrix(
                matrix_specifications, id),
            [](GlobalMatrix* const A)
            { NumLib::GlobalMatrixProvider::provider.releaseMatrix(*A); }};
}

//! Acquires a global vector from the provider, which takes it back when the
//! last owner releases it.
std::shared_ptr<GlobalVector> provideVector(
    MathLib::MatrixSpecifications const& matrix_specifications)
{
    std::size_t id = 0u;
    return {&NumLib::GlobalVectorProvider::provider.getVector(
                matrix_specifications, id),
            [](GlobalVector* const x)
            { NumLib::GlobalVectorProvider::provider.releaseVector(*x); }};
}
}  // namespace detail

namespace NumLib
{
TimeDiscretizedODESystem<ODESystemTag::FirstOrderImplicitQuasilinear,
                         NonlinearSolverTag::Newton>::
    TimeDiscretizedODESystem(const int process_id, ODE& ode,
                             TimeDisc& time_discretization)
    : _process_id(process_id),
      _ode(ode),
      _time_disc(time_discretization),
      _mat_trans(createMatrixTranslator<ODETag>(time_discretization)),
      _Jac(::detail::provideMatrix(_ode.getMatrixSpecifications(process_id))),
      _b(::detail::provideVector(_ode.getMatrixSpecifications(process_id)))
{
}

TimeDiscretizedODESystem<ODESystemTag::FirstOrderImplicitQuasilinear,
                         NonlinearSolverTag::Newton>::
    TimeDiscretizedODESystem(const int process_id, ODE& ode,
                             TimeDisc& time_discretization,
                             TimeDiscretizedODESystem const& storage_owner)
    : _process_id(process_id),
      _ode(ode),
      _time_disc(time_discretization),
      _mat_trans(createMatrixTranslator<ODETag>(time_discretization)),
      _Jac(storage_owner._Jac),
      _b(storage_owner._b)
{
    if (_ode.getMatrixSpecifications(_process_id) !=
        storage_owner._ode.getMatrixSpecifications(storage_owner._process_id))
    {
        OGS_FATAL(
            "Process {} cannot share the equation system storage of process "
            "{}: their matrix specifications differ.",
            _process_id, storage_owner._process_id);
    }
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
    try
    {
        _ode.preAssemble(t, dt, x_curr);
        _ode.assembleWithJacobian(t, dt, x_new_timestep, x_prev, process_id,
                                  *_b, *_Jac);
    }
    catch (AssemblyException const&)
    {
        MathLib::LinAlg::finalizeAssembly(*_b);
        MathLib::LinAlg::finalizeAssembly(*_Jac);
        throw;
    }
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

    MathLib::applyKnownSolution(
        Jac, res, minus_delta_x, ids, values,
        MathLib::DirichletBCApplicationMode::COMPLETE_MATRIX_UPDATE);
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
    MathLib::applyKnownSolution(
        Jac, res, x, ids, values,
        MathLib::DirichletBCApplicationMode::COMPLETE_MATRIX_UPDATE);
}

TimeDiscretizedODESystem<ODESystemTag::FirstOrderImplicitQuasilinear,
                         NonlinearSolverTag::Picard>::
    TimeDiscretizedODESystem(const int process_id, ODE& ode,
                             TimeDisc& time_discretization)
    : _process_id(process_id),
      _ode(ode),
      _time_disc(time_discretization),
      _mat_trans(createMatrixTranslator<ODETag>(time_discretization)),
      _M(::detail::provideMatrix(ode.getMatrixSpecifications(process_id))),
      _K(::detail::provideMatrix(ode.getMatrixSpecifications(process_id))),
      _b(::detail::provideVector(ode.getMatrixSpecifications(process_id)))
{
}

TimeDiscretizedODESystem<ODESystemTag::FirstOrderImplicitQuasilinear,
                         NonlinearSolverTag::Picard>::
    TimeDiscretizedODESystem(const int process_id, ODE& ode,
                             TimeDisc& time_discretization,
                             TimeDiscretizedODESystem const& storage_owner)
    : _process_id(process_id),
      _ode(ode),
      _time_disc(time_discretization),
      _mat_trans(createMatrixTranslator<ODETag>(time_discretization)),
      _M(storage_owner._M),
      _K(storage_owner._K),
      _b(storage_owner._b)
{
    if (_ode.getMatrixSpecifications(_process_id) !=
        storage_owner._ode.getMatrixSpecifications(storage_owner._process_id))
    {
        OGS_FATAL(
            "Process {} cannot share the equation system storage of process "
            "{}: their matrix specifications differ.",
            _process_id, storage_owner._process_id);
    }
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

    auto const finalize_assembly = [&]
    {
        LinAlg::finalizeAssembly(*_M);
        LinAlg::finalizeAssembly(*_K);
        LinAlg::finalizeAssembly(*_b);
    };

    try
    {
        _ode.preAssemble(t, dt, x_curr);
        _ode.assemble(t, dt, x_new_timestep, x_prev, process_id, *_M, *_K, *_b);
    }
    catch (AssemblyException const&)
    {
        finalize_assembly();
        throw;
    }
    finalize_assembly();
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
    applyKnownSolutionsPicard(
        GlobalMatrix& A, GlobalVector& rhs, GlobalVector& x,
        MathLib::DirichletBCApplicationMode const mode) const
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
    MathLib::applyKnownSolution(A, rhs, x, ids, values, mode);
}

}  // namespace NumLib
