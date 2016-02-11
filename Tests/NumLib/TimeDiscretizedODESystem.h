#pragma once

#include "ODESystem.h"
#include "NonlinearSystem.h"
#include "TimeDiscretization.h"
#include "MatrixTranslator.h"


template<NonlinearSolverTag NLTag_>
class TimeDiscretizedODESystemBase
        : public NonlinearSystem<NLTag_>
        , public InternalMatrixStorage
{
public:
    static constexpr NonlinearSolverTag NLTag  = NLTag_;

    virtual TimeDiscretization& getTimeDiscretization() = 0;
};


template<ODESystemTag ODETag, NonlinearSolverTag NLTag>
class TimeDiscretizedODESystem;


template<>
class TimeDiscretizedODESystem<ODESystemTag::FirstOrderImplicitQuasilinear,
                               NonlinearSolverTag::Newton> final
        : public TimeDiscretizedODESystemBase<NonlinearSolverTag::Newton>
{
public:
    static constexpr ODESystemTag ODETag = ODESystemTag::FirstOrderImplicitQuasilinear;

    using ODE = ODESystem<ODETag, NLTag>;
    using MatTrans = MatrixTranslator<ODETag>;


    explicit
    TimeDiscretizedODESystem(ODE& ode, TimeDiscretization& time_discretization, MatTrans& mat_trans)
        : _ode(ode)
        , _time_disc(time_discretization)
        , _mat_trans(mat_trans)
        , _Jac(ode.getMatrixSize(), ode.getMatrixSize())
        , _M(_Jac)
        , _K(_Jac)
        , _b(ode.getMatrixSize())
    {}

    /// begin INonlinearSystemNewton

    void assembleResidualNewton(const Vector &x_new_timestep) override
    {
        auto const  t      = _time_disc.getCurrentTime();
        auto const& x_curr = _time_disc.getCurrentX(x_new_timestep);

        _ode.assemble(t, x_curr, _M, _K, _b);
    }

    void assembleJacobian(const Vector &x_new_timestep) override
    {
        auto const  t        = _time_disc.getCurrentTime();
        auto const& x_curr   = _time_disc.getCurrentX(x_new_timestep);
        auto const  dxdot_dx = _time_disc.getCurrentXWeight();

        _ode.assembleJacobian(t, x_curr, dxdot_dx, _time_disc.getDxDx(), _Jac);
    }

    Vector getResidual(Vector const& x_new_timestep) override
    {
        return _mat_trans.getResidual(_M, _K, _b, x_new_timestep);
    }

    Matrix getJacobian() override
    {
        return _mat_trans.getJacobian(_Jac);
    }

    bool isLinear() const override
    {
        return _time_disc.isLinearTimeDisc() || _ode.isLinear();
    }

    TimeDiscretization& getTimeDiscretization() override {
        return _time_disc;
    }

    virtual void pushMatrices() const override
    {
        _mat_trans.pushMatrices(_M, _K, _b);
    }

private:
    ODE& _ode;
    TimeDiscretization& _time_disc;
    MatTrans& _mat_trans;

    Matrix _Jac;
    Matrix _M;
    Matrix _K;
    Vector _b;
};

template<>
class TimeDiscretizedODESystem<ODESystemTag::FirstOrderImplicitQuasilinear,
                               NonlinearSolverTag::Picard> final
        : public TimeDiscretizedODESystemBase<NonlinearSolverTag::Picard>
{
public:
    static constexpr ODESystemTag ODETag = ODESystemTag::FirstOrderImplicitQuasilinear;

    using ODE = ODESystem<ODETag, NLTag>;
    using MatTrans = MatrixTranslator<ODETag>;


    explicit
    TimeDiscretizedODESystem(ODE& ode, TimeDiscretization& time_discretization, MatTrans& mat_trans)
        : _ode(ode)
        , _time_disc(time_discretization)
        , _mat_trans(mat_trans)
        , _M(ode.getMatrixSize(), ode.getMatrixSize())
        , _K(_M)
        , _b(ode.getMatrixSize())
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

    TimeDiscretization& getTimeDiscretization() override {
        return _time_disc;
    }

    virtual void pushMatrices() const override
    {
        _mat_trans.pushMatrices(_M, _K, _b);
    }

private:
    ODE& _ode;
    TimeDiscretization& _time_disc;
    MatTrans& _mat_trans;

    Matrix _M;
    Matrix _K;
    Vector _b;
};
