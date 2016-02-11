#pragma once

#include "NonlinSolver.h"
#include "TimeDiscretization.h"

#include <memory>

template<typename Equation>
class MatrixTranslator;

template<>
class MatrixTranslator<IParabolicEquation>
{
public:
    virtual Matrix getA(Matrix const& M, Matrix const& K) const = 0;

    virtual Vector getRhs(const Matrix &M, const Matrix &K, const Vector& b) const = 0;

    virtual Vector getResidual(Matrix const& M, Matrix const& K, Vector const& b,
                               Vector const& x_new_timestep) const = 0;

    virtual Matrix getJacobian(Matrix Jac) const = 0;

    virtual void pushMatrices(Matrix const& M, Matrix const& K, Vector const& b)
    {
        (void) M; (void) K; (void) b;
    }

    virtual ~MatrixTranslator() = default;
};


template<typename Equation>
class MatrixTranslatorGeneral;

template<>
class MatrixTranslatorGeneral<IParabolicEquation>
        : public MatrixTranslator<IParabolicEquation>
{
public:
    MatrixTranslatorGeneral(ITimeDiscretization const& timeDisc)
        : _time_disc(timeDisc)
    {}

    Matrix getA(Matrix const& M, Matrix const& K) const override
    {
        auto const dxdot_dx = _time_disc.getCurrentXWeight();

        return M * dxdot_dx + K;
    }

    Vector getRhs(const Matrix &M, const Matrix &/*K*/, const Vector& b) const override
    {
        auto const& weighted_old_x = _time_disc.getWeightedOldX();

        return b + M * weighted_old_x;
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
    ITimeDiscretization const& _time_disc;
};


template<typename Equation>
class MatrixTranslatorForwardEuler;

template<>
class MatrixTranslatorForwardEuler<IParabolicEquation>
        : public MatrixTranslator<IParabolicEquation>
{
public:
    MatrixTranslatorForwardEuler(ForwardEuler const& timeDisc)
        : _fwd_euler(timeDisc)
    {}

    Matrix getA(Matrix const& M, Matrix const& /*K*/) const override
    {
        auto const dxdot_dx = _fwd_euler.getCurrentXWeight();

        return M * dxdot_dx;
    }

    Vector getRhs(const Matrix &M, const Matrix &K, const Vector& b) const override
    {
        auto const& weighted_old_x = _fwd_euler.getWeightedOldX();
        auto const& x_old          = _fwd_euler.getXOld();

        return b + M * weighted_old_x - K * x_old;
    }

    Vector getResidual(Matrix const& M, Matrix const& K, Vector const& b,
                       Vector const& x_new_timestep) const override
    {
        auto const  alpha  = _fwd_euler.getCurrentXWeight();
        auto const& x_curr = _fwd_euler.getCurrentX(x_new_timestep);
        auto const  x_old  = _fwd_euler.getWeightedOldX();
        auto const  x_dot  = alpha*x_new_timestep - x_old;

        return M * x_dot + K*x_curr - b;
    }

    Matrix getJacobian(Matrix Jac) const override
    {
        return Jac;
    }

private:
    ForwardEuler const& _fwd_euler;
};


template<typename Equation>
class MatrixTranslatorCrankNicolson;

template<>
class MatrixTranslatorCrankNicolson<IParabolicEquation>
        : public MatrixTranslator<IParabolicEquation>
{
public:
    MatrixTranslatorCrankNicolson(CrankNicolson const& timeDisc)
        : _crank_nicolson(timeDisc)
    {}

    Matrix getA(Matrix const& M, Matrix const& K) const override
    {
        auto const dxdot_dx = _crank_nicolson.getCurrentXWeight();
        auto const theta    = _crank_nicolson.getTheta();

        return theta * (M * dxdot_dx + K) + dxdot_dx * _M_bar;
    }

    Vector getRhs(const Matrix &M, const Matrix &/*K*/, const Vector& b) const override
    {
        auto const& weighted_old_x = _crank_nicolson.getWeightedOldX();
        auto const  theta          = _crank_nicolson.getTheta();

        return theta * (b + M * weighted_old_x) + _M_bar * weighted_old_x - _b_bar;
    }

    Vector getResidual(Matrix const& M, Matrix const& K, Vector const& b,
                       Vector const& x_new_timestep) const override
    {
        auto const  alpha  = _crank_nicolson.getCurrentXWeight();
        auto const& x_curr = _crank_nicolson.getCurrentX(x_new_timestep);
        auto const  x_old  = _crank_nicolson.getWeightedOldX();
        auto const  x_dot  = alpha*x_new_timestep - x_old;
        auto const  theta  = _crank_nicolson.getTheta();

        return theta * (M * x_dot + K*x_curr - b) + _M_bar * x_dot + _b_bar;
    }

    Matrix getJacobian(Matrix Jac) const override
    {
        auto const dxdot_dx = _crank_nicolson.getCurrentXWeight();
        auto const theta    = _crank_nicolson.getTheta();

        return theta * Jac + dxdot_dx * _M_bar;
    }

    void pushMatrices(Matrix const& M, Matrix const& K, Vector const& b) override
    {
        auto const theta = _crank_nicolson.getTheta();
        auto const x_old = _crank_nicolson.getXOld();

        _M_bar = (1.0-theta) * M;
        _b_bar = (1.0-theta) * (K * x_old - b);
    }

private:
    CrankNicolson const& _crank_nicolson;

    Matrix _M_bar;
    Vector _b_bar;
};

template<typename Equation>
std::unique_ptr<MatrixTranslator<Equation>>
createMatrixTranslator(ITimeDiscretization const& timeDisc)
{
    if (auto* fwd_euler = dynamic_cast<ForwardEuler const*>(&timeDisc)) {
        return std::unique_ptr<MatrixTranslator<Equation>>(
                new MatrixTranslatorForwardEuler<Equation>(*fwd_euler));
    } else if (auto* crank = dynamic_cast<CrankNicolson const*>(&timeDisc)) {
        return std::unique_ptr<MatrixTranslator<Equation>>(
                new MatrixTranslatorCrankNicolson<Equation>(*crank));
    } else {
        return std::unique_ptr<MatrixTranslator<Equation>>(
                new MatrixTranslatorGeneral<Equation>(timeDisc));
    }
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
                             MatrixTranslator<IParabolicEquation>& mat_trans)
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

    /// end INonlinearSystemNewton

    ITimeDiscretization& getTimeDiscretization() {
        return _time_disc;
    }

private:
    // from IParabolicEquation
    virtual void pushMatrices() const override
    {
        _mat_trans.pushMatrices(_M, _K, _b);
    }


    IFirstOrderImplicitOde<NonlinearSolverTag::Newton>& _ode;
    ITimeDiscretization& _time_disc;
    MatrixTranslator<IParabolicEquation>& _mat_trans;

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
                             MatrixTranslator<IParabolicEquation>& mat_trans)
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
        return _mat_trans.getA(_M, _K);
    }

    Vector getRhs() override
    {
        return _mat_trans.getRhs(_M, _K, _b);
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
    virtual void pushMatrices() const override
    {
        _mat_trans.pushMatrices(_M, _K, _b);
    }


    IFirstOrderImplicitOde<NonlinearSolverTag::Picard>& _ode;
    ITimeDiscretization& _time_disc;
    MatrixTranslator<IParabolicEquation>& _mat_trans;

    Matrix _M;
    Matrix _K;
    Vector _b;
};
