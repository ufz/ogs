#pragma once

#include "ODESystem.h"
#include "NonlinearSystem.h"
#include "TimeDiscretization.h"
#include "MatrixTranslator.h"


namespace NumLib
{

template<typename Matrix, typename Vector, NonlinearSolverTag NLTag_>
class TimeDiscretizedODESystemBase
        : public NonlinearSystem<Matrix, Vector, NLTag_>
        , public InternalMatrixStorage
{
public:
    virtual TimeDiscretization<Vector>& getTimeDiscretization() = 0;
};


template<typename Matrix, typename Vector, ODESystemTag ODETag, NonlinearSolverTag NLTag>
class TimeDiscretizedODESystem;


template<typename Matrix, typename Vector>
class TimeDiscretizedODESystem<Matrix, Vector,
                               ODESystemTag::FirstOrderImplicitQuasilinear,
                               NonlinearSolverTag::Newton> final
        : public TimeDiscretizedODESystemBase<Matrix, Vector, NonlinearSolverTag::Newton>
{
public:
    static const ODESystemTag ODETag = ODESystemTag::FirstOrderImplicitQuasilinear;

    using ODE = ODESystem<Matrix, Vector, ODETag, NonlinearSolverTag::Newton>;
    using MatTrans = MatrixTranslator<Matrix, Vector, ODETag>;
    using TimeDisc = TimeDiscretization<Vector>;


    explicit
    TimeDiscretizedODESystem(ODE& ode, TimeDisc& time_discretization, MatTrans& mat_trans)
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
        namespace BLAS = MathLib::BLAS;

        auto const  t        = _time_disc.getCurrentTime();
        auto const& x_curr   = _time_disc.getCurrentX(x_new_timestep);
        auto const  dxdot_dx = _time_disc.getCurrentXWeight();

        _time_disc.getWeightedOldX(_xdot);
        BLAS::axpby(_xdot, dxdot_dx, -1.0, x_new_timestep);

        _ode.assembleJacobian(t, x_curr, _xdot, dxdot_dx, _time_disc.getDxDx(), _Jac);
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
    ODE& _ode;
    TimeDisc& _time_disc;
    MatTrans& _mat_trans;

    Matrix _Jac;
    Matrix _M;
    Matrix _K;
    Vector _b;
    Vector _xdot; // cache only
};

template<typename Matrix, typename Vector>
class TimeDiscretizedODESystem<Matrix, Vector,
                               ODESystemTag::FirstOrderImplicitQuasilinear,
                               NonlinearSolverTag::Picard> final
        : public TimeDiscretizedODESystemBase<Matrix, Vector, NonlinearSolverTag::Picard>
{
public:
    static const ODESystemTag ODETag = ODESystemTag::FirstOrderImplicitQuasilinear;

    using ODE = ODESystem<Matrix, Vector, ODETag, NonlinearSolverTag::Picard>;
    using MatTrans = MatrixTranslator<Matrix, Vector, ODETag>;
    using TimeDisc = TimeDiscretization<Vector>;


    explicit
    TimeDiscretizedODESystem(ODE& ode, TimeDisc& time_discretization, MatTrans& mat_trans)
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

    TimeDisc& getTimeDiscretization() override {
        return _time_disc;
    }

    virtual void pushMatrices() const override
    {
        _mat_trans.pushMatrices(_M, _K, _b);
    }

private:
    ODE& _ode;
    TimeDisc& _time_disc;
    MatTrans& _mat_trans;

    Matrix _M;
    Matrix _K;
    Vector _b;
};

}
