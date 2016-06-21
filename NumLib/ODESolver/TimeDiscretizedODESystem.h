/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef NUMLIB_TIMEDISCRETIZEDODESYSTEM_H
#define NUMLIB_TIMEDISCRETIZEDODESYSTEM_H

#include <memory>

#include "MathLib/LinAlg/ApplyKnownSolution.h"
#include "MathLib/LinAlg/UnifiedMatrixSetters.h"
#include "NumLib/IndexValueVector.h"

#include "MatrixTranslator.h"
#include "NonlinearSystem.h"
#include "ODESystem.h"
#include "TimeDiscretization.h"

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
}
}

namespace NumLib
{
//! \addtogroup ODESolver
//! @{

/*! A NonlinearSystem together with some TimeDiscretization scheme.
 *
 * This is the interface of an ODE towards the TimeLoop.
 * This interface is abstract, it represents any type of first-order ODE.
 *
 * \tparam Matrix the type of matrices occuring in the linearization of the ODE.
 * \tparam Vector the type of the solution vector of the ODE.
 * \tparam NLTag  a tag indicating the method used for resolving nonlinearities.
 */
template <NonlinearSolverTag NLTag>
class TimeDiscretizedODESystemBase : public NonlinearSystem<NLTag>,
                                     public InternalMatrixStorage
{
public:
    //! Exposes the used time discretization scheme.
    virtual TimeDiscretization& getTimeDiscretization() = 0;
};

/*! A NonlinearSystem together with some TimeDiscretization scheme.
 *
 * This class represents a specific type of first-order ODE, as indicated by
 * \c ODETag.
 *
 * \tparam Matrix the type of matrices occuring in the linearization of the ODE.
 * \tparam Vector the type of the solution vector of the ODE.
 * \tparam ODETag a tag indicating the type of ODE.
 * \tparam NLTag  a tag indicating the method used for resolving nonlinearities.
 */
template <ODESystemTag ODETag, NonlinearSolverTag NLTag>
class TimeDiscretizedODESystem;

/*! Time discretized first order implicit quasi-linear ODE;
 *  to be solved using the Newton-Raphson method for resolving nonlinearities.
 *
 * \tparam Matrix the type of matrices occuring in the linearization of the ODE.
 * \tparam Vector the type of the solution vector of the ODE.
 *
 * \see ODESystemTag::FirstOrderImplicitQuasilinear
 */
template <>
class TimeDiscretizedODESystem<ODESystemTag::FirstOrderImplicitQuasilinear,
                               NonlinearSolverTag::Newton>
    final : public TimeDiscretizedODESystemBase<NonlinearSolverTag::Newton>
{
public:
    //! A tag indicating the type of ODE.
    static const ODESystemTag ODETag =
        ODESystemTag::FirstOrderImplicitQuasilinear;

    //! The type of ODE.
    using ODE = ODESystem<ODETag, NonlinearSolverTag::Newton>;
    //! The auxiliary class that computes the matrix/vector used by the
    //! nonlinear solver.
    using MatTrans = MatrixTranslator<ODETag>;
    //! A shortcut for a general time discretization scheme
    using TimeDisc = TimeDiscretization;

    /*! Constructs a new instance.
     *
     * \param ode the ODE to be wrapped.
     * \param time_discretization the time discretization to be used.
     */
    explicit TimeDiscretizedODESystem(ODE& ode, TimeDisc& time_discretization)
        : _ode(ode),
          _time_disc(time_discretization),
          _mat_trans(createMatrixTranslator<ODETag>(time_discretization))
    {
        _Jac = &MathLib::GlobalMatrixProvider<GlobalMatrix>::provider.getMatrix(
            _ode.getMatrixSpecifications(), _Jac_id);
        _M = &MathLib::GlobalMatrixProvider<GlobalMatrix>::provider.getMatrix(
            _ode.getMatrixSpecifications(), _M_id);
        _K = &MathLib::GlobalMatrixProvider<GlobalMatrix>::provider.getMatrix(
            _ode.getMatrixSpecifications(), _K_id);
        _b = &MathLib::GlobalVectorProvider<GlobalVector>::provider.getVector(
            _ode.getMatrixSpecifications(), _b_id);
    }

    ~TimeDiscretizedODESystem()
    {
        MathLib::GlobalMatrixProvider<GlobalMatrix>::provider.releaseMatrix(
            *_Jac);
        MathLib::GlobalMatrixProvider<GlobalMatrix>::provider.releaseMatrix(
            *_M);
        MathLib::GlobalMatrixProvider<GlobalMatrix>::provider.releaseMatrix(
            *_K);
        MathLib::GlobalVectorProvider<GlobalVector>::provider.releaseVector(
            *_b);
    }

    void assembleResidualNewton(const GlobalVector& x_new_timestep) override
    {
        namespace BLAS = MathLib::BLAS;

        auto const t = _time_disc.getCurrentTime();
        auto const& x_curr = _time_disc.getCurrentX(x_new_timestep);

        _M->setZero();
        _K->setZero();
        _b->setZero();

        _ode.assemble(t, x_curr, *_M, *_K, *_b);

        BLAS::finalizeAssembly(*_M);
        BLAS::finalizeAssembly(*_K);
        BLAS::finalizeAssembly(*_b);
    }

    void assembleJacobian(const GlobalVector& x_new_timestep) override
    {
        namespace BLAS = MathLib::BLAS;

        auto const t = _time_disc.getCurrentTime();
        auto const& x_curr = _time_disc.getCurrentX(x_new_timestep);
        auto const dxdot_dx = _time_disc.getNewXWeight();
        auto const dx_dx = _time_disc.getDxDx();

        auto& xdot =
            MathLib::GlobalVectorProvider<GlobalVector>::provider.getVector(
                _xdot_id);
        _time_disc.getXdot(x_new_timestep, xdot);

        _Jac->setZero();

        _ode.assembleJacobian(t, x_curr, xdot, dxdot_dx, *_M, dx_dx, *_K,
                              *_Jac);

        MathLib::BLAS::finalizeAssembly(*_Jac);

        MathLib::GlobalVectorProvider<GlobalVector>::provider.releaseVector(
            xdot);
    }

    void getResidual(GlobalVector const& x_new_timestep,
                     GlobalVector& res) const override
    {
        // TODO Maybe the duplicate calculation of xdot here and in
        //      assembleJacobian can be optimuized. However, that would make
        //      the interface a bit more fragile.
        auto& xdot =
            MathLib::GlobalVectorProvider<GlobalVector>::provider.getVector(
                _xdot_id);
        _time_disc.getXdot(x_new_timestep, xdot);

        _mat_trans->computeResidual(*_M, *_K, *_b, x_new_timestep, xdot, res);

        MathLib::GlobalVectorProvider<GlobalVector>::provider.releaseVector(
            xdot);
    }

    void getJacobian(GlobalMatrix& Jac) const override
    {
        _mat_trans->computeJacobian(*_Jac, Jac);
    }

    void applyKnownSolutions(GlobalVector& x) const override
    {
        ::detail::applyKnownSolutions(
            _ode.getKnownSolutions(_time_disc.getCurrentTime()), x);
    }

    void applyKnownSolutionsNewton(GlobalMatrix& Jac, GlobalVector& res,
                                   GlobalVector& minus_delta_x) override
    {
        auto const* known_solutions =
            _ode.getKnownSolutions(_time_disc.getCurrentTime());

        if (known_solutions)
        {
            std::vector<double> values;

            for (auto const& bc : *known_solutions)
            {
                // TODO this is the quick and dirty and bad performance
                // solution.
                values.resize(bc.values.size(), 0.0);

                // TODO maybe it would be faster to apply all at once
                MathLib::applyKnownSolution(Jac, res, minus_delta_x, bc.ids,
                                            values);
            }
        }
    }

    bool isLinear() const override
    {
        return _time_disc.isLinearTimeDisc() || _ode.isLinear();
    }

    void preIteration(const unsigned iter, GlobalVector const& x) override
    {
        _ode.preIteration(iter, x);
    }

    IterationResult postIteration(GlobalVector const& x) override
    {
        return _ode.postIteration(x);
    }

    void pushMatrices() const override
    {
        _mat_trans->pushMatrices(*_M, *_K, *_b);
    }

    TimeDisc& getTimeDiscretization() override { return _time_disc; }
    MathLib::MatrixSpecifications getMatrixSpecifications() const override
    {
        return _ode.getMatrixSpecifications();
    }

private:
    ODE& _ode;             //!< ode the ODE being wrapped
    TimeDisc& _time_disc;  //!< the time discretization to being used

    //! the object used to compute the matrix/vector for the nonlinear solver
    std::unique_ptr<MatTrans> _mat_trans;

    GlobalMatrix* _Jac;  //!< the Jacobian of the residual
    GlobalMatrix* _M;    //!< Matrix \f$ M \f$.
    GlobalMatrix* _K;    //!< Matrix \f$ K \f$.
    GlobalVector* _b;    //!< Matrix \f$ b \f$.

    std::size_t _Jac_id = 0u;  //!< ID of the \c _Jac matrix.
    std::size_t _M_id = 0u;    //!< ID of the \c _M matrix.
    std::size_t _K_id = 0u;    //!< ID of the \c _K matrix.
    std::size_t _b_id = 0u;    //!< ID of the \c _b vector.

    //! ID of the vector storing xdot in intermediate computations.
    mutable std::size_t _xdot_id = 0u;
};

/*! Time discretized first order implicit quasi-linear ODE;
 *  to be solved using the Picard fixpoint iteration method for resolving
 * nonlinearities.
 *
 * \tparam Matrix the type of matrices occuring in the linearization of the ODE.
 * \tparam Vector the type of the solution vector of the ODE.
 *
 * \see ODESystemTag::FirstOrderImplicitQuasilinear
 */
template <>
class TimeDiscretizedODESystem<ODESystemTag::FirstOrderImplicitQuasilinear,
                               NonlinearSolverTag::Picard>
    final : public TimeDiscretizedODESystemBase<NonlinearSolverTag::Picard>
{
public:
    //! A tag indicating the type of ODE.
    static const ODESystemTag ODETag =
        ODESystemTag::FirstOrderImplicitQuasilinear;

    //! The type of ODE.
    using ODE = ODESystem<ODETag, NonlinearSolverTag::Picard>;
    //! The auxiliary class that computes the matrix/vector used by the
    //! nonlinear solver.
    using MatTrans = MatrixTranslator<ODETag>;
    //! A shortcut for a general time discretization scheme
    using TimeDisc = TimeDiscretization;

    /*! Constructs a new instance.
     *
     * \param ode the ODE to be wrapped.
     * \param time_discretization the time discretization to be used.
     */
    explicit TimeDiscretizedODESystem(ODE& ode, TimeDisc& time_discretization)
        : _ode(ode),
          _time_disc(time_discretization),
          _mat_trans(createMatrixTranslator<ODETag>(time_discretization))
    {
        _M = &MathLib::GlobalMatrixProvider<GlobalMatrix>::provider.getMatrix(
            ode.getMatrixSpecifications(), _M_id);
        _K = &MathLib::GlobalMatrixProvider<GlobalMatrix>::provider.getMatrix(
            ode.getMatrixSpecifications(), _K_id);
        _b = &MathLib::GlobalVectorProvider<GlobalVector>::provider.getVector(
            ode.getMatrixSpecifications(), _b_id);
    }

    ~TimeDiscretizedODESystem()
    {
        MathLib::GlobalMatrixProvider<GlobalMatrix>::provider.releaseMatrix(
            *_M);
        MathLib::GlobalMatrixProvider<GlobalMatrix>::provider.releaseMatrix(
            *_K);
        MathLib::GlobalVectorProvider<GlobalVector>::provider.releaseVector(
            *_b);
    }

    void assembleMatricesPicard(const GlobalVector& x_new_timestep) override
    {
        namespace BLAS = MathLib::BLAS;

        auto const t = _time_disc.getCurrentTime();
        auto const& x_curr = _time_disc.getCurrentX(x_new_timestep);

        _M->setZero();
        _K->setZero();
        _b->setZero();

        _ode.assemble(t, x_curr, *_M, *_K, *_b);

        BLAS::finalizeAssembly(*_M);
        BLAS::finalizeAssembly(*_K);
        BLAS::finalizeAssembly(*_b);
    }

    void getA(GlobalMatrix& A) const override
    {
        _mat_trans->computeA(*_M, *_K, A);
    }

    void getRhs(GlobalVector& rhs) const override
    {
        _mat_trans->computeRhs(*_M, *_K, *_b, rhs);
    }

    void applyKnownSolutions(GlobalVector& x) const override
    {
        ::detail::applyKnownSolutions(
            _ode.getKnownSolutions(_time_disc.getCurrentTime()), x);
    }

    void applyKnownSolutionsPicard(GlobalMatrix& A, GlobalVector& rhs,
                                   GlobalVector& x) override
    {
        auto const* known_solutions =
            _ode.getKnownSolutions(_time_disc.getCurrentTime());

        if (known_solutions)
        {
            using IndexType =
                typename MathLib::MatrixVectorTraits<GlobalMatrix>::Index;
            std::vector<IndexType> ids;
            std::vector<double> values;
            for (auto const& bc : *known_solutions)
            {
                std::copy(bc.ids.cbegin(), bc.ids.cend(),
                          std::back_inserter(ids));
                std::copy(bc.values.cbegin(), bc.values.cend(),
                          std::back_inserter(values));
            }
            MathLib::applyKnownSolution(A, rhs, x, ids, values);
        }
    }

    bool isLinear() const override
    {
        return _time_disc.isLinearTimeDisc() || _ode.isLinear();
    }

    void preIteration(const unsigned iter, GlobalVector const& x) override
    {
        _ode.preIteration(iter, x);
    }

    IterationResult postIteration(GlobalVector const& x) override
    {
        return _ode.postIteration(x);
    }

    void pushMatrices() const override
    {
        _mat_trans->pushMatrices(*_M, *_K, *_b);
    }

    TimeDisc& getTimeDiscretization() override { return _time_disc; }
    MathLib::MatrixSpecifications getMatrixSpecifications() const override
    {
        return _ode.getMatrixSpecifications();
    }

private:
    ODE& _ode;             //!< ode the ODE being wrapped
    TimeDisc& _time_disc;  //!< the time discretization to being used

    //! the object used to compute the matrix/vector for the nonlinear solver
    std::unique_ptr<MatTrans> _mat_trans;

    GlobalMatrix* _M;  //!< Matrix \f$ M \f$.
    GlobalMatrix* _K;  //!< Matrix \f$ K \f$.
    GlobalVector* _b;  //!< Matrix \f$ b \f$.

    std::size_t _M_id = 0u;  //!< ID of the \c _M matrix.
    std::size_t _K_id = 0u;  //!< ID of the \c _K matrix.
    std::size_t _b_id = 0u;  //!< ID of the \c _b vector.
};

//! @}
}

#endif  // NUMLIB_TIMEDISCRETIZEDODESYSTEM_H
