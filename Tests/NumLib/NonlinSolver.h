#pragma once

#include <Eigen/SparseCore>

using Matrix = Eigen::SparseMatrix<double>;
using Vector = Eigen::VectorXd;


class INonlinearSystemNewton
{
public:
    virtual void assembleResidualNewton(Vector const& x) = 0;
    virtual void assembleJacobian(Vector const& x) = 0;
    virtual Vector getResidual(Vector const& x) = 0;
    virtual Matrix getJacobian() = 0;
    virtual ~INonlinearSystemNewton() = default;
};

class INonlinearSystemPicard
{
public:
    virtual void assembleMatricesPicard(Vector const& x) = 0;
    virtual Matrix getA() = 0;
    virtual Vector getRhs() = 0;
    virtual ~INonlinearSystemPicard() = default;
};



class NonlinearSolverNewton final
{
public:
    NonlinearSolverNewton(double const tol, const unsigned maxiter)
        : _tol(tol)
        , _maxiter(maxiter)
    {}

    void solve(INonlinearSystemNewton& sys, Vector& x);

private:
    const double _tol;
    const unsigned _maxiter;

    Vector _minus_delta_x;
};


class NonlinearSolverPicard final
{
public:
    NonlinearSolverPicard(double const tol, const unsigned maxiter)
        : _tol(tol)
        , _maxiter(maxiter)
    {}

    void solve(INonlinearSystemPicard& sys, Vector& x);

private:
    const double _tol;
    const unsigned _maxiter;

    Vector _x_new;
};



class ITimeDiscretization
{
public:
    virtual void pushState(const double t, Vector const& x) = 0;
    virtual void pushMatrices() = 0;
    virtual void setCurrentTime(const double t, const double delta_t) = 0;
    virtual double getCurrentTime() const = 0;

    // \dot x === alpha * x - x_old
    virtual double getCurrentXWeight() = 0; // = alpha
    virtual Vector getWeightedOldX() = 0; // = x_old

    ~ITimeDiscretization() = default;
};



class IFirstOrderImplicitOde
{
public:
    virtual void assemble(const double t, Vector const& x,
                          Matrix& M, Matrix& K, Vector& b) = 0;


    virtual ~IFirstOrderImplicitOde() = default;
};

class IFirstOrderImplicitOdeNewton : public IFirstOrderImplicitOde
{
public:
    // TODO: \dot x contribution
    virtual void assembleJacobian(const double t, Vector const& x,
                                  const double dxdot_dx,
                                  Matrix& Jac) = 0;
};



enum class NonlinearSolverTag : bool { Picard, Newton };


template<NonlinearSolverTag NonlinearSolver, typename TimeDisc>
struct TimeDiscretizedODESystem;

template<typename TimeDisc>
class TimeDiscretizedODESystem<NonlinearSolverTag::Newton, TimeDisc> final
        : public TimeDisc
        , public INonlinearSystemNewton
{
public:
    TimeDiscretizedODESystem(IFirstOrderImplicitOdeNewton& ode)
        : _ode(ode)
    {}

    void assembleResidualNewton(const Vector &x) override
    {
        auto const t = TimeDisc::getCurrentTime();
        _ode.assemble(t, x, _M, _K, _b);
    }

    void assembleJacobian(const Vector &x) override
    {
        auto const t = TimeDisc::getCurrentTime();
        auto const dxdot_dx = TimeDisc::getCurrentXWeight();
        _ode.assembleJacobian(t, x, dxdot_dx, _Jac);
    }

    Vector getResidual(Vector const& x)
    {
        auto const alpha = TimeDisc::getCurrentXWeight();
        auto const x_old = TimeDisc::getWeightedOldX();
        return _M * (alpha*x - x_old) + _K*x - _b;
    }

    Matrix getJacobian() {
        return _Jac;
    }

private:
    IFirstOrderImplicitOdeNewton& _ode;

    Matrix _Jac;
    Matrix _M;
    Matrix _K;
    Vector _b;
};

template<typename TimeDisc>
class TimeDiscretizedODESystem<NonlinearSolverTag::Picard, TimeDisc> final
        : public TimeDisc
        , public INonlinearSystemPicard
{
public:
    TimeDiscretizedODESystem(IFirstOrderImplicitOde& ode)
        : _ode(ode)
    {}

    void assembleMatricesPicard(const Vector &x) override
    {
        auto const t = TimeDisc::getCurrentTime();
        _ode.assemble(t, x, _M, _K, _b);
    }

    Matrix getA() override {
        return _M * TimeDisc::getCurrentXWeight() + _K;
    }

    Vector getRhs() override {
        return _b + _M * TimeDisc::getWeightedOldX();
    }

private:
    IFirstOrderImplicitOde& _ode;

    Matrix _M;
    Matrix _K;
    Vector _b;
};

//

