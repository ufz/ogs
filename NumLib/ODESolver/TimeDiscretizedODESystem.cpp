/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "TimeDiscretizedODESystem.h"

#include "MathLib/LinAlg/ApplyKnownSolution.h"
#include "MathLib/LinAlg/UnifiedMatrixSetters.h"
#include "NumLib/IndexValueVector.h"

namespace detail
{
//! Applies known solutions to the solution vector \c x.
template <typename Solutions, typename Vector>
void applyKnownSolutions(std::vector<Solutions> const* const known_solutions,
                         Vector& x)
{
    if (!known_solutions)
        return;

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
}

namespace NumLib
{
TimeDiscretizedODESystem<ODESystemTag::FirstOrderImplicitQuasilinear,
                         NonlinearSolverTag::Newton>::
    TimeDiscretizedODESystem(ODE& ode, TimeDisc& time_discretization)
    : _ode(ode),
      _time_disc(time_discretization),
      _mat_trans(createMatrixTranslator<ODETag>(time_discretization))
{
    _Jac = &NumLib::GlobalMatrixProvider::provider.getMatrix(
        _ode.getMatrixSpecifications(), _Jac_id);
    _M = &NumLib::GlobalMatrixProvider::provider.getMatrix(
        _ode.getMatrixSpecifications(), _M_id);
    _K = &NumLib::GlobalMatrixProvider::provider.getMatrix(
        _ode.getMatrixSpecifications(), _K_id);
    _b = &NumLib::GlobalVectorProvider::provider.getVector(
        _ode.getMatrixSpecifications(), _b_id);
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
    assemble(const GlobalVector& x_new_timestep)
{
    namespace LinAlg = MathLib::LinAlg;

    auto const t = _time_disc.getCurrentTime();
    auto const& x_curr = _time_disc.getCurrentX(x_new_timestep);
    auto const dxdot_dx = _time_disc.getNewXWeight();
    auto const dx_dx = _time_disc.getDxDx();

    auto& xdot = NumLib::GlobalVectorProvider::provider.getVector(_xdot_id);
    _time_disc.getXdot(x_new_timestep, xdot);

    _M->setZero();
    _K->setZero();
    _b->setZero();
    _Jac->setZero();

    _ode.assembleWithJacobian(t, x_curr, xdot, dxdot_dx, dx_dx, *_M, *_K, *_b,
                              *_Jac);

    LinAlg::finalizeAssembly(*_M);
    LinAlg::finalizeAssembly(*_K);
    LinAlg::finalizeAssembly(*_b);
    MathLib::LinAlg::finalizeAssembly(*_Jac);

    NumLib::GlobalVectorProvider::provider.releaseVector(xdot);
}

void TimeDiscretizedODESystem<
    ODESystemTag::FirstOrderImplicitQuasilinear,
    NonlinearSolverTag::Newton>::getResidual(GlobalVector const& x_new_timestep,
                                             GlobalVector& res) const
{
    // TODO Maybe the duplicate calculation of xdot here and in assembleJacobian
    //      can be optimuized. However, that would make the interface a bit more
    //      fragile.
    auto& xdot = NumLib::GlobalVectorProvider::provider.getVector(_xdot_id);
    _time_disc.getXdot(x_new_timestep, xdot);

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
    NonlinearSolverTag::Newton>::applyKnownSolutions(GlobalVector& x) const
{
    ::detail::applyKnownSolutions(
        _ode.getKnownSolutions(_time_disc.getCurrentTime()), x);
}

void TimeDiscretizedODESystem<ODESystemTag::FirstOrderImplicitQuasilinear,
                              NonlinearSolverTag::Newton>::
    applyKnownSolutionsNewton(GlobalMatrix& Jac, GlobalVector& res,
                              GlobalVector& minus_delta_x)
{
    auto const* known_solutions =
        _ode.getKnownSolutions(_time_disc.getCurrentTime());

    if (!known_solutions || known_solutions->empty())
        return;

    using IndexType = MathLib::MatrixVectorTraits<GlobalMatrix>::Index;
    std::vector<IndexType> ids;
    for (auto const& bc : *known_solutions)
    {
        std::copy(bc.ids.cbegin(), bc.ids.cend(), std::back_inserter(ids));
    }

    // For the Newton method the values must be zero
    std::vector<double> values(ids.size(), 0);
    MathLib::applyKnownSolution(Jac, res, minus_delta_x, ids, values);
}

TimeDiscretizedODESystem<ODESystemTag::FirstOrderImplicitQuasilinear,
                         NonlinearSolverTag::Picard>::
    TimeDiscretizedODESystem(ODE& ode, TimeDisc& time_discretization)
    : _ode(ode),
      _time_disc(time_discretization),
      _mat_trans(createMatrixTranslator<ODETag>(time_discretization))
{
    _M = &NumLib::GlobalMatrixProvider::provider.getMatrix(
        ode.getMatrixSpecifications(), _M_id);
    _K = &NumLib::GlobalMatrixProvider::provider.getMatrix(
        ode.getMatrixSpecifications(), _K_id);
    _b = &NumLib::GlobalVectorProvider::provider.getVector(
        ode.getMatrixSpecifications(), _b_id);
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
    assemble(const GlobalVector& x_new_timestep)
{
    namespace LinAlg = MathLib::LinAlg;

    auto const t = _time_disc.getCurrentTime();
    auto const& x_curr = _time_disc.getCurrentX(x_new_timestep);

    _M->setZero();
    _K->setZero();
    _b->setZero();

    _ode.assemble(t, x_curr, *_M, *_K, *_b);

    LinAlg::finalizeAssembly(*_M);
    LinAlg::finalizeAssembly(*_K);
    LinAlg::finalizeAssembly(*_b);
}

void TimeDiscretizedODESystem<
    ODESystemTag::FirstOrderImplicitQuasilinear,
    NonlinearSolverTag::Picard>::applyKnownSolutions(GlobalVector& x) const
{
    ::detail::applyKnownSolutions(
        _ode.getKnownSolutions(_time_disc.getCurrentTime()), x);
}

void TimeDiscretizedODESystem<
    ODESystemTag::FirstOrderImplicitQuasilinear,
    NonlinearSolverTag::Picard>::applyKnownSolutionsPicard(GlobalMatrix& A,
                                                           GlobalVector& rhs,
                                                           GlobalVector& x)
{
    auto const* known_solutions =
        _ode.getKnownSolutions(_time_disc.getCurrentTime());

    if (known_solutions)
    {
        using IndexType = MathLib::MatrixVectorTraits<GlobalMatrix>::Index;
        std::vector<IndexType> ids;
        std::vector<double> values;
        for (auto const& bc : *known_solutions)
        {
            std::copy(bc.ids.cbegin(), bc.ids.cend(), std::back_inserter(ids));
            std::copy(bc.values.cbegin(), bc.values.cend(),
                      std::back_inserter(values));
        }
        MathLib::applyKnownSolution(A, rhs, x, ids, values);
    }
}

}  // NumLib
