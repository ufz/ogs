#pragma once

#include <memory>

#include "ODETypes.h"
#include "TimeDiscretization.h"

#include "BLAS.h"


template<ODESystemTag ODETag>
class MatrixTranslator;

template<>
class MatrixTranslator<ODESystemTag::FirstOrderImplicitQuasilinear>
{
public:
    virtual Matrix getA(Matrix const& M, Matrix const& K) const = 0;

    virtual Vector getRhs(const Matrix &M, const Matrix &K, const Vector& b) const = 0;

    virtual Vector getResidual(Matrix const& M, Matrix const& K, Vector const& b,
                               Vector const& x_new_timestep) const = 0;

    virtual Matrix getJacobian(Matrix const& Jac) const = 0;

    // needed for Crank-Nicolson
    virtual void pushMatrices(Matrix const& M, Matrix const& K, Vector const& b)
    {
        (void) M; (void) K; (void) b;
    }

    virtual ~MatrixTranslator() = default;
};


template<ODESystemTag ODETag>
class MatrixTranslatorGeneral;

template<>
class MatrixTranslatorGeneral<ODESystemTag::FirstOrderImplicitQuasilinear>
        : public MatrixTranslator<ODESystemTag::FirstOrderImplicitQuasilinear>
{
public:
    MatrixTranslatorGeneral(TimeDiscretization const& timeDisc)
        : _time_disc(timeDisc)
    {}

    Matrix getA(Matrix const& M, Matrix const& K) const override
    {
        auto const dxdot_dx = _time_disc.getCurrentXWeight();

        // A = M * dxdot_dx + K
        Matrix A(M);
        BLAS::aypx(A, dxdot_dx, K);

        return A;
    }

    Vector getRhs(const Matrix &M, const Matrix &/*K*/, const Vector& b) const override
    {
        auto const& weighted_old_x = _time_disc.getWeightedOldX();

        // rhs = M * weighted_old_x + b
        Vector rhs;
        BLAS::matMultAdd(M, weighted_old_x, b, rhs);

        return rhs;
    }

    Vector getResidual(Matrix const& M, Matrix const& K, Vector const& b,
                       Vector const& x_new_timestep) const override
    {
        auto const  alpha  = _time_disc.getCurrentXWeight();
        auto const& x_curr = _time_disc.getCurrentX(x_new_timestep);
        auto const  x_old  = _time_disc.getWeightedOldX();
        auto const  x_dot  = alpha*x_new_timestep - x_old;

        // res = M * x_dot + K * x_curr - b
        Vector res;
        BLAS::matMult(M, x_dot, res);
        BLAS::matMultAdd(K, x_curr, res, res);
        BLAS::axpy(res, -1.0, b);

        return res;
    }

    Matrix getJacobian(Matrix const& Jac) const override
    {
        return Jac;
    }

private:
    TimeDiscretization const& _time_disc;
};


template<ODESystemTag ODETag>
class MatrixTranslatorForwardEuler;

template<>
class MatrixTranslatorForwardEuler<ODESystemTag::FirstOrderImplicitQuasilinear>
        : public MatrixTranslator<ODESystemTag::FirstOrderImplicitQuasilinear>
{
public:
    MatrixTranslatorForwardEuler(ForwardEuler const& timeDisc)
        : _fwd_euler(timeDisc)
    {}

    Matrix getA(Matrix const& M, Matrix const& /*K*/) const override
    {
        auto const dxdot_dx = _fwd_euler.getCurrentXWeight();

        // A = M * dxdot_dx
        Matrix A(M);
        BLAS::scale(A, dxdot_dx);

        return A;
    }

    Vector getRhs(const Matrix &M, const Matrix &K, const Vector& b) const override
    {
        auto const& weighted_old_x = _fwd_euler.getWeightedOldX();
        auto const& x_old          = _fwd_euler.getXOld();

        // rhs = b + M * weighted_old_x - K * x_old
        Vector rhs;
        BLAS::matMult(K, x_old, rhs); // rhs = K * x_old
        BLAS::aypx(rhs, -1.0, b);     // rhs = b - K * x_old
        BLAS::matMultAdd(M, weighted_old_x, rhs, rhs); // rhs += M * weighted_old_x

        return rhs;
    }

    Vector getResidual(Matrix const& M, Matrix const& K, Vector const& b,
                       Vector const& x_new_timestep) const override
    {
        auto const  alpha  = _fwd_euler.getCurrentXWeight();
        auto const& x_curr = _fwd_euler.getCurrentX(x_new_timestep);
        auto const  x_old  = _fwd_euler.getWeightedOldX();
        auto const  x_dot  = alpha*x_new_timestep - x_old;

        // res = M * x_dot + K * x_curr - b
        Vector res;
        BLAS::matMult(M, x_dot, res);
        BLAS::matMultAdd(K, x_curr, res, res);
        BLAS::axpy(res, -1.0, b);

        return res;
    }

    Matrix getJacobian(Matrix const& Jac) const override
    {
        return Jac;
    }

private:
    ForwardEuler const& _fwd_euler;
};


template<ODESystemTag ODETag>
class MatrixTranslatorCrankNicolson;

template<>
class MatrixTranslatorCrankNicolson<ODESystemTag::FirstOrderImplicitQuasilinear>
        : public MatrixTranslator<ODESystemTag::FirstOrderImplicitQuasilinear>
{
public:
    MatrixTranslatorCrankNicolson(CrankNicolson const& timeDisc)
        : _crank_nicolson(timeDisc)
    {}

    Matrix getA(Matrix const& M, Matrix const& K) const override
    {
        auto const dxdot_dx = _crank_nicolson.getCurrentXWeight();
        auto const theta    = _crank_nicolson.getTheta();

        // A = theta * (M * dxdot_dx + K) + dxdot_dx * _M_bar
        Matrix A(M);
        BLAS::aypx(A, dxdot_dx, K); // A = M * dxdot_dx + K

        BLAS::scale(A, theta); // A *= theta
        BLAS::axpy(A, dxdot_dx, _M_bar); // A += dxdot_dx * _M_bar

        return A;
    }

    Vector getRhs(const Matrix &M, const Matrix &/*K*/, const Vector& b) const override
    {
        auto const& weighted_old_x = _crank_nicolson.getWeightedOldX();
        auto const  theta          = _crank_nicolson.getTheta();

        // rhs = theta * (b + M * weighted_old_x) + _M_bar * weighted_old_x - _b_bar;
        Vector rhs;
        BLAS::matMultAdd(M, weighted_old_x, b, rhs); // rhs = b + M * weighted_old_x

        BLAS::scale(rhs, theta); // rhs *= theta
        BLAS::matMultAdd(_M_bar, weighted_old_x, rhs, rhs); // rhs += _M_bar * weighted_old_x
        BLAS::axpy(rhs, -1.0, _b_bar); // rhs -= b

        return rhs;
    }

    Vector getResidual(Matrix const& M, Matrix const& K, Vector const& b,
                       Vector const& x_new_timestep) const override
    {
        auto const  alpha  = _crank_nicolson.getCurrentXWeight();
        auto const& x_curr = _crank_nicolson.getCurrentX(x_new_timestep);
        auto const  x_old  = _crank_nicolson.getWeightedOldX();
        auto const  x_dot  = alpha*x_new_timestep - x_old;
        auto const  theta  = _crank_nicolson.getTheta();

        // res = theta * (M * x_dot + K*x_curr - b) + _M_bar * x_dot + _b_bar
        Vector res;
        BLAS::matMult(M, x_dot, res); // res = M * x_dot
        BLAS::matMultAdd(K, x_curr, res, res); // res += K * x_curr
        BLAS::axpy(res, -1.0, b); // res = M * x_dot + K * x_curr - b

        BLAS::aypx(res, theta, _b_bar); // res = res * theta + _b_bar
        BLAS::matMultAdd(_M_bar, x_dot, res, res); // rs += _M_bar * x_dot

        return res;
    }

    Matrix getJacobian(Matrix const& Jac) const override
    {
        auto const dxdot_dx = _crank_nicolson.getCurrentXWeight();
        auto const theta    = _crank_nicolson.getTheta();

        // J = theta * Jac + dxdot_dx * _M_bar
        Matrix J(Jac);
        BLAS::scale(J, theta);
        BLAS::axpy(J, dxdot_dx, _M_bar);

        return J;
    }

    void pushMatrices(Matrix const& M, Matrix const& K, Vector const& b) override
    {
        auto const theta = _crank_nicolson.getTheta();
        auto const x_old = _crank_nicolson.getXOld();

        // _M_bar = (1.0-theta) * M;
        BLAS::copy(M, _M_bar);
        BLAS::scale(_M_bar, 1.0-theta);

        // _b_bar = (1.0-theta) * (K * x_old - b)
        BLAS::matMult(K, x_old, _b_bar);
        BLAS::axpy(_b_bar, -1.0, b);
        BLAS::scale(_b_bar, 1.0-theta);
    }

private:
    CrankNicolson const& _crank_nicolson;

    Matrix _M_bar;
    Vector _b_bar;
};


template<ODESystemTag ODETag>
std::unique_ptr<MatrixTranslator<ODETag>>
createMatrixTranslator(TimeDiscretization const& timeDisc)
{
    if (auto* fwd_euler = dynamic_cast<ForwardEuler const*>(&timeDisc)) {
        return std::unique_ptr<MatrixTranslator<ODETag>>(
                new MatrixTranslatorForwardEuler<ODETag>(*fwd_euler));
    } else if (auto* crank = dynamic_cast<CrankNicolson const*>(&timeDisc)) {
        return std::unique_ptr<MatrixTranslator<ODETag>>(
                new MatrixTranslatorCrankNicolson<ODETag>(*crank));
    } else {
        return std::unique_ptr<MatrixTranslator<ODETag>>(
                new MatrixTranslatorGeneral<ODETag>(timeDisc));
    }
}
