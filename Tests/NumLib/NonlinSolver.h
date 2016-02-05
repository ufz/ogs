#pragma once

#include <Eigen/SparseCore>
#include <logog/include/logog.hpp>

using Matrix = Eigen::SparseMatrix<double>;
using Vector = Eigen::VectorXd;


enum class NonlinearSolverTag : bool { Picard, Newton };


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


template<NonlinearSolverTag NLTag>
class NonlinearSolver;

template<>
class NonlinearSolver<NonlinearSolverTag::Newton> final
{
public:
    NonlinearSolver(double const tol, const unsigned maxiter)
        : _tol(tol)
        , _maxiter(maxiter)
    {}

    void solve(INonlinearSystemNewton& sys, Vector& x);

private:
    const double _tol;
    const unsigned _maxiter;

    Vector _minus_delta_x;
};

template<>
class NonlinearSolver<NonlinearSolverTag::Picard> final
{
public:
    NonlinearSolver(double const tol, const unsigned maxiter)
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


template<NonlinearSolverTag NLTag>
class IFirstOrderImplicitOde;

template<>
class IFirstOrderImplicitOde<NonlinearSolverTag::Picard>
{
public:
    virtual void assemble(const double t, Vector const& x,
                          Matrix& M, Matrix& K, Vector& b) = 0;


    virtual ~IFirstOrderImplicitOde() = default;
};

template<>
class IFirstOrderImplicitOde<NonlinearSolverTag::Newton>
        : public IFirstOrderImplicitOde<NonlinearSolverTag::Picard>
{
public:
    // TODO: \dot x contribution
    virtual void assembleJacobian(const double t, Vector const& x,
                                  const double dxdot_dx,
                                  Matrix& Jac) = 0;
};




template<NonlinearSolverTag NLTag, typename TimeDisc>
struct TimeDiscretizedODESystem;

template<typename TimeDisc>
class TimeDiscretizedODESystem<NonlinearSolverTag::Newton, TimeDisc> final
        : public TimeDisc
        , public INonlinearSystemNewton
{
public:
    TimeDiscretizedODESystem(IFirstOrderImplicitOde<NonlinearSolverTag::Newton>& ode)
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
    IFirstOrderImplicitOde<NonlinearSolverTag::Newton>& _ode;

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
    TimeDiscretizedODESystem(IFirstOrderImplicitOde<NonlinearSolverTag::Picard>& ode)
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
    IFirstOrderImplicitOde<NonlinearSolverTag::Picard>& _ode;

    Matrix _M;
    Matrix _K;
    Vector _b;
};


template<NonlinearSolverTag NLTag, typename TimeDisc>
class TimeLoop
{
public:
    using TDiscODESys = TimeDiscretizedODESystem<NLTag, TimeDisc>;
    using NLSolver = NonlinearSolver<NLTag>;

    TimeLoop(TDiscODESys& ode_sys, NLSolver nonlinear_solver)
        : _ode_sys(ode_sys)
        , _nonlinear_solver(nonlinear_solver)
    {}

    void loop(const double t0, const Vector x0, const double t_end, const double delta_t);

private:
    TDiscODESys& _ode_sys;
    NLSolver& _nonlinear_solver;
};


template<NonlinearSolverTag NLTag, typename TimeDisc>
void
TimeLoop<NLTag, TimeDisc>::
loop(const double t0, const Vector x0, const double t_end, const double delta_t)
{
    Vector x(x0); // solution vector

    _ode_sys.pushState(t0, x0); // push IC

    for (double t=t0+delta_t; t<=t_end; t+=delta_t)
    {
        INFO("time: %e, delta_t: %e", t, delta_t);
        _ode_sys.setCurrentTime(t, delta_t);

        _nonlinear_solver.solve(_ode_sys, x);

        _ode_sys.pushState(t, x);
        _ode_sys.pushMatrices();

        INFO("x[0] = %10g, x[1] = %10g", x[0], x[1]);
        // auto const v = ode1_solution(t);
        // INFO("x[0] = %10g, x[1] = %10g (analytical solution)", v[0], v[1]);
    }
}



//

