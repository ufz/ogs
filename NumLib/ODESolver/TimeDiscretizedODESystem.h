/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <memory>

#include "MatrixTranslator.h"
#include "NonlinearSystem.h"
#include "ODESystem.h"
#include "TimeDiscretization.h"

namespace NumLib
{
//! \addtogroup ODESolver
//! @{

/*! A NonlinearSystem together with some TimeDiscretization scheme.
 *
 * This is the interface of an ODE towards the TimeLoop.
 * This interface is abstract, it represents any type of first-order ODE.
 *
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
 * \tparam ODETag a tag indicating the type of ODE.
 * \tparam NLTag  a tag indicating the method used for resolving nonlinearities.
 */
template <ODESystemTag ODETag, NonlinearSolverTag NLTag>
class TimeDiscretizedODESystem;

/*! Time discretized first order implicit quasi-linear ODE;
 *  to be solved using the Newton-Raphson method for resolving nonlinearities.
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
     * \param process_id ID of the ODE to be solved.
     * \param ode the ODE to be wrapped.
     * \param time_discretization the time discretization to be used.
     */
    explicit TimeDiscretizedODESystem(const int process_id, ODE& ode,
                                      TimeDisc& time_discretization);

    ~TimeDiscretizedODESystem() override;

    void assemble(const GlobalVector& x_new_timestep) override;

    void getResidual(GlobalVector const& x_new_timestep,
                     GlobalVector& res) const override;

    void getJacobian(GlobalMatrix& Jac) const override;

    void applyKnownSolutions(GlobalVector& x) const override;

    void applyKnownSolutionsNewton(GlobalMatrix& Jac, GlobalVector& res,
                                   GlobalVector& minus_delta_x,
                                   GlobalVector& x) override;

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
    MathLib::MatrixSpecifications getMatrixSpecifications(
        const int process_id) const override
    {
        return _ode.getMatrixSpecifications(process_id);
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

    //! Constructs a new instance.
    explicit TimeDiscretizedODESystem(const int process_id, ODE& ode,
                                      TimeDisc& time_discretization);

    ~TimeDiscretizedODESystem() override;

    void assemble(const GlobalVector& x_new_timestep) override;

    void getA(GlobalMatrix& A) const override
    {
        _mat_trans->computeA(*_M, *_K, A);
    }

    void getRhs(GlobalVector& rhs) const override
    {
        _mat_trans->computeRhs(*_M, *_K, *_b, rhs);
    }

    void applyKnownSolutions(GlobalVector& x) const override;

    void applyKnownSolutionsPicard(GlobalMatrix& A, GlobalVector& rhs,
                                   GlobalVector& x) override;

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
    MathLib::MatrixSpecifications getMatrixSpecifications(
        const int process_id) const override
    {
        return _ode.getMatrixSpecifications(process_id);
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
