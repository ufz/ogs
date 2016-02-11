#pragma once

#include "NonlinSolver.h"
#include "TimeDiscretization.h"

#include <memory>

template<typename Equation>
class MatrixTranslatorBase;

template<>
class MatrixTranslatorBase<IParabolicEquation>
{
public:
    virtual Matrix getA(Matrix const& M, Matrix const& K) const = 0;

    virtual Vector getRhs(const Matrix &M, const Matrix &K, const Vector& b) const = 0;

    virtual Vector getResidual(Matrix const& M, Matrix const& K, Vector const& b,
                               Vector const& x_new_timestep) const = 0;

    virtual Matrix getJacobian(Matrix Jac) const = 0;

    virtual ~MatrixTranslatorBase() = default;
};


template<typename TimeDisc, typename Equation>
class MatrixTranslator;

// TODO: don't need timedisc template param
template<typename TimeDisc>
class MatrixTranslator<TimeDisc, IParabolicEquation>
        : public MatrixTranslatorBase<IParabolicEquation>
{
public:
    MatrixTranslator(TimeDisc const& timeDisc)
        : _time_disc(timeDisc)
    {}

    Matrix getA(Matrix const& M, Matrix const& K) const override
    {
        auto const dxdot_dx = _time_disc.getCurrentXWeight();
        auto const dx_dx    = _time_disc.getDxDx();
        return M * dxdot_dx + K * dx_dx;
    }

    Vector getRhs(const Matrix &M, const Matrix &K, const Vector& b) const override
    {
        auto const& weighted_old_x = _time_disc.getWeightedOldX();
        auto const& x_old = b; // TODO: fix
        auto const  dx_dx = _time_disc.getDxDx();
        return b + M * weighted_old_x - (1.0-dx_dx) * K * x_old;
    }

    Vector getResidual(Matrix const& M, Matrix const& K, Vector const& b,
                       Vector const& x_new_timestep) const override
    {
        auto const  alpha  = _time_disc.getCurrentXWeight();
        auto const& x_curr = _time_disc.getCurrentX(x_new_timestep);
        auto const  x_old  = _time_disc.getWeightedOldX();
        auto const  x_dot  = alpha*x_new_timestep - x_old;

        return M * x_dot + K*x_curr - b;
    }

    Matrix getJacobian(Matrix Jac) const override
    {
        return Jac;
    }

private:
    TimeDisc const& _time_disc;
};

template<typename Equation>
std::unique_ptr<MatrixTranslatorBase<Equation>>
createMatrixTranslator(ITimeDiscretization const& timeDisc)
{
    return std::unique_ptr<MatrixTranslatorBase<Equation>>(
            new MatrixTranslator<ITimeDiscretization, Equation>(timeDisc));
}



template<NonlinearSolverTag NLTag>
struct TimeDiscretizedODESystem;

template<>
class TimeDiscretizedODESystem<NonlinearSolverTag::Newton> final
        : public INonlinearSystemNewton
        , public IParabolicEquation
{
public:
    explicit
    TimeDiscretizedODESystem(IFirstOrderImplicitOde<NonlinearSolverTag::Newton>& ode,
                             ITimeDiscretization& time_discretization,
                             MatrixTranslatorBase<IParabolicEquation>& mat_trans)
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
        _time_disc.adjustMatrix(_Jac);
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

    /// end INonlinearSystemNewton

    ITimeDiscretization& getTimeDiscretization() {
        return _time_disc;
    }

private:
    // from IParabolicEquation
    virtual void getMatrices(Matrix const*& M, Matrix const*& K,
                             Vector const*& b) const override
    {
        M = &_M;
        K = &_K;
        b = &_b;
    }


    IFirstOrderImplicitOde<NonlinearSolverTag::Newton>& _ode;
    ITimeDiscretization& _time_disc;
    MatrixTranslatorBase<IParabolicEquation>& _mat_trans;

    Matrix _Jac;
    Matrix _M;
    Matrix _K;
    Vector _b;
};

template<>
class TimeDiscretizedODESystem<NonlinearSolverTag::Picard> final
        : public INonlinearSystemPicard
        , public IParabolicEquation
{
public:
    explicit
    TimeDiscretizedODESystem(IFirstOrderImplicitOde<NonlinearSolverTag::Picard>& ode,
                             ITimeDiscretization& time_discretization,
                             MatrixTranslatorBase<IParabolicEquation>& mat_trans)
        : _ode(ode)
        , _time_disc(time_discretization)
        , _mat_trans(mat_trans)
        , _M(ode.getMatrixSize(), ode.getMatrixSize())
        , _K(_M)
        , _b(ode.getMatrixSize())
    {}

    /// begin INonlinearSystemNewton

    void assembleMatricesPicard(const Vector &x_new_timestep) override
    {
        auto const  t      = _time_disc.getCurrentTime();
        auto const& x_curr = _time_disc.getCurrentX(x_new_timestep);

        _ode.assemble(t, x_curr, _M, _K, _b);
    }

    Matrix getA() override
    {
        return _time_disc.getA(_M, _K);
    }

    Vector getRhs() override
    {
        return _time_disc.getRhs(_M, _K, _b);
    }

    bool isLinear() const override
    {
        return _time_disc.isLinearTimeDisc() || _ode.isLinear();
    }

    /// end INonlinearSystemPicard

    ITimeDiscretization& getTimeDiscretization() {
        return _time_disc;
    }

private:
    // from IParabolicEquation
    virtual void getMatrices(Matrix const*& M, Matrix const*& K,
                             Vector const*& b) const override
    {
        M = &_M;
        K = &_K;
        b = &_b;
    }


    IFirstOrderImplicitOde<NonlinearSolverTag::Picard>& _ode;
    ITimeDiscretization& _time_disc;
    MatrixTranslatorBase<IParabolicEquation>& _mat_trans;

    Matrix _M;
    Matrix _K;
    Vector _b;
};
