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

#include "ODESystem.h"
#include "NonlinearSystem.h"
#include "TimeDiscretization.h"
#include "MatrixTranslator.h"


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
template<typename Matrix, typename Vector, NonlinearSolverTag NLTag>
class TimeDiscretizedODESystemBase
        : public NonlinearSystem<Matrix, Vector, NLTag>
        , public InternalMatrixStorage
{
public:
    virtual TimeDiscretization<Vector>& getTimeDiscretization() = 0;
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
template<typename Matrix, typename Vector, ODESystemTag ODETag, NonlinearSolverTag NLTag>
class TimeDiscretizedODESystem;


/*! Time discretized first order implicit quasi-linear ODE;
 *  to be solved using the Newton-Raphson method for resolving nonlinearities.
 *
 * \tparam Matrix the type of matrices occuring in the linearization of the ODE.
 * \tparam Vector the type of the solution vector of the ODE.
 *
 * \see ODESystemTag::FirstOrderImplicitQuasilinear
 */
template<typename Matrix, typename Vector>
class TimeDiscretizedODESystem<Matrix, Vector,
                               ODESystemTag::FirstOrderImplicitQuasilinear,
                               NonlinearSolverTag::Newton> final
        : public TimeDiscretizedODESystemBase<Matrix, Vector, NonlinearSolverTag::Newton>
{
public:
    //! A tag indicating the type of ODE.
    static const ODESystemTag ODETag = ODESystemTag::FirstOrderImplicitQuasilinear;

    //! The type of ODE.
    using ODE = ODESystem<Matrix, Vector, ODETag, NonlinearSolverTag::Newton>;
    //! The auxiliary class that computes the matrix/vector used by the nonlinear solver.
    using MatTrans = MatrixTranslator<Matrix, Vector, ODETag>;
    //! A shortcut for a general time discretization scheme
    using TimeDisc = TimeDiscretization<Vector>;

    /*! Constructs a new instance.
     *
     * \param ode the ODE to be wrapped.
     * \param time_discretization the time discretization to be used.
     * \param mat_trans the object used to compute the matrix/vector for the nonlinear solver.
     */
    explicit
    TimeDiscretizedODESystem(ODE& ode, TimeDisc& time_discretization, MatTrans& mat_trans)
        : _ode(ode)
        , _time_disc(time_discretization)
        , _mat_trans(mat_trans)
        , _Jac(ode.getNumEquations(), ode.getNumEquations())
        , _M(_Jac)
        , _K(_Jac)
        , _b(ode.getNumEquations())
    {}

    void assembleResidualNewton(const Vector &x_new_timestep) override
    {
        auto const  t      = _time_disc.getCurrentTime();
        auto const& x_curr = _time_disc.getCurrentX(x_new_timestep);

        _ode.assemble(t, x_curr, _M, _K, _b);
    }

    void assembleJacobian(const Vector &x_new_timestep) override
    {
        namespace BLAS = MathLib::BLAS;

        auto const  t        = _time_disc.getCurrentTime();
        auto const& x_curr   = _time_disc.getCurrentX(x_new_timestep);
        auto const  dxdot_dx = _time_disc.getCurrentXWeight();
        auto const  dx_dx    = _time_disc.getDxDx();
        _time_disc.getXdot(x_new_timestep, _xdot);

        _ode.assembleJacobian(t, x_curr, _xdot, dxdot_dx, dx_dx, _Jac);
    }

    void getResidual(Vector const& x_new_timestep, Vector& res) override
    {
        _mat_trans.getResidual(_M, _K, _b, x_new_timestep, res);
    }

    void getJacobian(Matrix& Jac) override
    {
        _mat_trans.getJacobian(_Jac, Jac);
    }

    bool isLinear() const override
    {
        return _time_disc.isLinearTimeDisc() || _ode.isLinear();
    }

    TimeDisc& getTimeDiscretization() override {
        return _time_disc;
    }

    virtual void pushMatrices() const override
    {
        _mat_trans.pushMatrices(_M, _K, _b);
    }

private:
    ODE& _ode;            //!< ode the ODE being wrapped
    TimeDisc& _time_disc; //!< the time discretization to being used
    MatTrans& _mat_trans; //!< the object used to compute the matrix/vector for the nonlinear solver

    Matrix _Jac;  //!< the Jacobian of the residual
    Matrix _M;    //!< Matrix \f$ M \f$.
    Matrix _K;    //!< Matrix \f$ K \f$.
    Vector _b;    //!< Matrix \f$ b \f$.

    Vector _xdot; //!< Used to cache \f$ \dot x \f$. \todo Save some memory.
};


/*! Time discretized first order implicit quasi-linear ODE;
 *  to be solved using the Picard fixpoint iteration method for resolving nonlinearities.
 *
 * \tparam Matrix the type of matrices occuring in the linearization of the ODE.
 * \tparam Vector the type of the solution vector of the ODE.
 *
 * \see ODESystemTag::FirstOrderImplicitQuasilinear
 */
template<typename Matrix, typename Vector>
class TimeDiscretizedODESystem<Matrix, Vector,
                               ODESystemTag::FirstOrderImplicitQuasilinear,
                               NonlinearSolverTag::Picard> final
        : public TimeDiscretizedODESystemBase<Matrix, Vector, NonlinearSolverTag::Picard>
{
public:
    //! A tag indicating the type of ODE.
    static const ODESystemTag ODETag = ODESystemTag::FirstOrderImplicitQuasilinear;

    //! The type of ODE.
    using ODE = ODESystem<Matrix, Vector, ODETag, NonlinearSolverTag::Picard>;
    //! The auxiliary class that computes the matrix/vector used by the nonlinear solver.
    using MatTrans = MatrixTranslator<Matrix, Vector, ODETag>;
    //! A shortcut for a general time discretization scheme
    using TimeDisc = TimeDiscretization<Vector>;

    /*! Constructs a new instance.
     *
     * \param ode the ODE to be wrapped.
     * \param time_discretization the time discretization to be used.
     * \param mat_trans the object used to compute the matrix/vector for the nonlinear solver.
     */
    explicit
    TimeDiscretizedODESystem(ODE& ode, TimeDisc& time_discretization, MatTrans& mat_trans)
        : _ode(ode)
        , _time_disc(time_discretization)
        , _mat_trans(mat_trans)
        , _M(ode.getNumEquations(), ode.getNumEquations())
        , _K(_M)
        , _b(ode.getNumEquations())
    {}

    void assembleMatricesPicard(const Vector &x_new_timestep) override
    {
        auto const  t      = _time_disc.getCurrentTime();
        auto const& x_curr = _time_disc.getCurrentX(x_new_timestep);

        _ode.assemble(t, x_curr, _M, _K, _b);
    }

    void getA(Matrix& A) override
    {
        _mat_trans.getA(_M, _K, A);
    }

    void getRhs(Vector& rhs) override
    {
        _mat_trans.getRhs(_M, _K, _b, rhs);
    }

    bool isLinear() const override
    {
        return _time_disc.isLinearTimeDisc() || _ode.isLinear();
    }

    TimeDisc& getTimeDiscretization() override {
        return _time_disc;
    }

    virtual void pushMatrices() const override
    {
        _mat_trans.pushMatrices(_M, _K, _b);
    }

private:
    ODE& _ode;            //!< ode the ODE being wrapped
    TimeDisc& _time_disc; //!< the time discretization to being used
    MatTrans& _mat_trans; //!< the object used to compute the matrix/vector for the nonlinear solver

    Matrix _M;    //!< Matrix \f$ M \f$.
    Matrix _K;    //!< Matrix \f$ K \f$.
    Vector _b;    //!< Matrix \f$ b \f$.
};

//! @}

}

#endif // NUMLIB_TIMEDISCRETIZEDODESYSTEM_H
