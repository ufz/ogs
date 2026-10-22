// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "TimeDiscretizedODESystem.h"

#include <range/v3/view/zip.hpp>

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

//! \return The differences x_i - g_i between \c x and the known solutions'
//! values g at the known solutions' ids.
template <typename Index, typename Vector>
std::vector<double> differencesToKnownSolutions(
    std::vector<Index> const& ids,
    std::vector<double> const& values,
    Vector const& x)
{
    std::vector<double> differences;
    differences.reserve(ids.size());
    for (auto const& [id, value] : ranges::views::zip(ids, values))
    {
        differences.push_back(x[id] - value);
    }
    return differences;
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

    // Concatenate the ids and values; clear() keeps the capacity.
    _known_solutions_ids.clear();
    _known_solutions_values.clear();
    if (_known_solutions)
    {
        for (auto const& bc : *_known_solutions)
        {
            _known_solutions_ids.insert(_known_solutions_ids.end(),
                                        bc.ids.begin(), bc.ids.end());
            _known_solutions_values.insert(_known_solutions_values.end(),
                                           bc.values.begin(), bc.values.end());
        }
    }
}

void TimeDiscretizedODESystem<
    ODESystemTag::FirstOrderImplicitQuasilinear,
    NonlinearSolverTag::Newton>::applyKnownSolutions(GlobalVector& x) const
{
    ::detail::applyKnownSolutions(_known_solutions, x);
}

double TimeDiscretizedODESystem<ODESystemTag::FirstOrderImplicitQuasilinear,
                                NonlinearSolverTag::Newton>::
    getDirichletFreeResidual(GlobalVector const& x, GlobalVector const& res,
                             GlobalVector& w)
{
    auto& differences = NumLib::GlobalVectorProvider::provider.getVector(
        res, _known_solutions_differences_id);

    // V_k = x_k - g_k at the Dirichlet DOFs, zero elsewhere.
    differences.setZero();
    MathLib::setVector(differences, _known_solutions_ids,
                       ::detail::differencesToKnownSolutions(
                           _known_solutions_ids, _known_solutions_values, x));
    MathLib::LinAlg::finalizeAssembly(differences);
    double const max_deviation =
        MathLib::LinAlg::norm(differences, MathLib::VecNormType::INFINITY_N);

    // w = res - J V
    if (max_deviation > 0)
    {
        MathLib::LinAlg::matMult(*_Jac, differences, w);
        MathLib::LinAlg::aypx(w, -1.0, res);
    }
    else
    {
        MathLib::LinAlg::copy(res, w);
    }
    NumLib::GlobalVectorProvider::provider.releaseVector(differences);

    // Restrict to the free DOFs.
    MathLib::setVector(w, _known_solutions_ids,
                       std::vector<double>(_known_solutions_ids.size(), 0.0));
    MathLib::LinAlg::finalizeAssembly(w);

    return max_deviation;
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

    // minus_delta_x will be set to the difference between the current value
    // and the Dirichlet BC value, s.t. the Newton step imposes the BC values.
    // After a full step it vanishes up to round-off.
    auto const values = ::detail::differencesToKnownSolutions(
        _known_solutions_ids, _known_solutions_values, x);

    MathLib::applyKnownSolution(
        Jac, res, minus_delta_x, _known_solutions_ids, values,
        MathLib::DirichletBCApplicationMode::COMPLETE_MATRIX_UPDATE);
}

void TimeDiscretizedODESystem<ODESystemTag::FirstOrderImplicitQuasilinear,
                              NonlinearSolverTag::Newton>::
    applyKnownSolutionsPETScSNES(GlobalMatrix& Jac, GlobalVector& res,
                                 GlobalVector const& x,
                                 GlobalVector& scratch) const
{
    if (!_known_solutions)
    {
        return;
    }

    // res_i = x_i - g_i is zero for an iterate fulfilling the Dirichlet BCs.
    // The iterate itself must not be passed to applyKnownSolution(), which
    // overwrites its third argument with the given values.
    auto const values = ::detail::differencesToKnownSolutions(
        _known_solutions_ids, _known_solutions_values, x);
    MathLib::applyKnownSolution(
        Jac, res, scratch, _known_solutions_ids, values,
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
