#pragma once

#include <memory>

#include "Types.h"
#include "TimeDiscretization.h"

#include "MathLib/LinAlg/BLAS.h"

namespace NumLib
{

//! \addtogroup ODESolver
//! @{

template<typename Matrix, typename Vector, ODESystemTag ODETag>
class MatrixTranslator;

template<typename Matrix, typename Vector>
class MatrixTranslator<Matrix, Vector, ODESystemTag::FirstOrderImplicitQuasilinear>
{
public:
    virtual void getA(Matrix const& M, Matrix const& K, Matrix& A) const = 0;

    virtual void getRhs(const Matrix &M, const Matrix &K, const Vector& b, Vector& rhs) const = 0;

    virtual void getResidual(Matrix const& M, Matrix const& K, Vector const& b,
                             Vector const& x_new_timestep,
                             Vector& res) const = 0;

    virtual void getJacobian(Matrix const& Jac_in, Matrix& Jac_out) const = 0;

    // needed for Crank-Nicolson
    virtual void pushMatrices(Matrix const& M, Matrix const& K, Vector const& b)
    {
        (void) M; (void) K; (void) b;
    }

    virtual ~MatrixTranslator() = default;
};


template<typename Matrix, typename Vector, ODESystemTag ODETag>
class MatrixTranslatorGeneral;

template<typename Matrix, typename Vector>
class MatrixTranslatorGeneral<Matrix, Vector, ODESystemTag::FirstOrderImplicitQuasilinear>
        : public MatrixTranslator<Matrix, Vector, ODESystemTag::FirstOrderImplicitQuasilinear>
{
public:
    MatrixTranslatorGeneral(TimeDiscretization<Vector> const& timeDisc)
        : _time_disc(timeDisc)
    {}

    void getA(Matrix const& M, Matrix const& K, Matrix& A) const override
    {
        namespace BLAS = MathLib::BLAS;

        auto const dxdot_dx = _time_disc.getCurrentXWeight();

        // A = M * dxdot_dx + K
        BLAS::copy(M, A);
        BLAS::aypx(A, dxdot_dx, K);
    }

    void getRhs(const Matrix &M, const Matrix &/*K*/, const Vector& b, Vector& rhs) const override
    {
        namespace BLAS = MathLib::BLAS;

        _time_disc.getWeightedOldX(_tmp);

        // rhs = M * weighted_old_x + b
        BLAS::matMultAdd(M, _tmp, b, rhs);
    }

    void getResidual(Matrix const& M, Matrix const& K, Vector const& b,
                     Vector const& x_new_timestep, Vector& res) const override
    {
        namespace BLAS = MathLib::BLAS;

        auto const  alpha  = _time_disc.getCurrentXWeight();
        auto const& x_curr = _time_disc.getCurrentX(x_new_timestep);

        // x_dot  = alpha*x_new_timestep - x_old
        _time_disc.getWeightedOldX(_tmp);
        BLAS::axpby(_tmp, alpha, -1.0, x_new_timestep);

        // res = M * x_dot + K * x_curr - b
        BLAS::matMult(M, _tmp, res); // the local vector x_dot seems to be necessary because of this multiplication
        BLAS::matMultAdd(K, x_curr, res, res);
        BLAS::axpy(res, -1.0, b);
    }

    void getJacobian(Matrix const& Jac_in, Matrix& Jac_out) const override
    {
        namespace BLAS = MathLib::BLAS;

        BLAS::copy(Jac_in, Jac_out);
    }

private:
    TimeDiscretization<Vector> const& _time_disc;
    mutable Vector _tmp; // used to store intermediate calculation results
};


template<typename Matrix, typename Vector, ODESystemTag ODETag>
class MatrixTranslatorForwardEuler;

template<typename Matrix, typename Vector>
class MatrixTranslatorForwardEuler<Matrix, Vector, ODESystemTag::FirstOrderImplicitQuasilinear>
        : public MatrixTranslator<Matrix, Vector, ODESystemTag::FirstOrderImplicitQuasilinear>
{
public:
    MatrixTranslatorForwardEuler(ForwardEuler<Vector> const& timeDisc)
        : _fwd_euler(timeDisc)
    {}

    void getA(Matrix const& M, Matrix const& /*K*/, Matrix& A) const override
    {
        namespace BLAS = MathLib::BLAS;

        auto const dxdot_dx = _fwd_euler.getCurrentXWeight();

        // A = M * dxdot_dx
        BLAS::copy(M, A);
        BLAS::scale(A, dxdot_dx);
    }

    void getRhs(const Matrix &M, const Matrix &K, const Vector& b, Vector& rhs) const override
    {
        namespace BLAS = MathLib::BLAS;

        _fwd_euler.getWeightedOldX(_tmp);

        auto const& x_old          = _fwd_euler.getXOld();

        // rhs = b + M * weighted_old_x - K * x_old
        BLAS::matMult(K, x_old, rhs); // rhs = K * x_old
        BLAS::aypx(rhs, -1.0, b);     // rhs = b - K * x_old
        BLAS::matMultAdd(M, _tmp, rhs, rhs); // rhs += M * weighted_old_x
    }

    void getResidual(Matrix const& M, Matrix const& K, Vector const& b,
                     Vector const& x_new_timestep,
                     Vector& res) const override
    {
        namespace BLAS = MathLib::BLAS;

        auto const  alpha  = _fwd_euler.getCurrentXWeight();
        auto const& x_curr = _fwd_euler.getCurrentX(x_new_timestep);

        // x_dot  = alpha*x_new_timestep - x_old
        _fwd_euler.getWeightedOldX(_tmp);
        BLAS::axpby(_tmp, alpha, -1.0, x_new_timestep);

        // res = M * x_dot + K * x_curr - b
        BLAS::matMult(M, _tmp, res); // the local vector x_dot seems to be necessary because of this multiplication
        BLAS::matMultAdd(K, x_curr, res, res);
        BLAS::axpy(res, -1.0, b);
    }

    void getJacobian(Matrix const& Jac_in, Matrix& Jac_out) const override
    {
        namespace BLAS = MathLib::BLAS;

        BLAS::copy(Jac_in, Jac_out);
    }

private:
    ForwardEuler<Vector> const& _fwd_euler;
    mutable Vector _tmp; // used to store intermediate calculation results
};


template<typename Matrix, typename Vector, ODESystemTag ODETag>
class MatrixTranslatorCrankNicolson;

template<typename Matrix, typename Vector>
class MatrixTranslatorCrankNicolson<Matrix, Vector, ODESystemTag::FirstOrderImplicitQuasilinear>
        : public MatrixTranslator<Matrix, Vector, ODESystemTag::FirstOrderImplicitQuasilinear>
{
public:
    MatrixTranslatorCrankNicolson(CrankNicolson<Vector> const& timeDisc)
        : _crank_nicolson(timeDisc)
    {}

    void getA(Matrix const& M, Matrix const& K, Matrix& A) const override
    {
        namespace BLAS = MathLib::BLAS;

        auto const dxdot_dx = _crank_nicolson.getCurrentXWeight();
        auto const theta    = _crank_nicolson.getTheta();

        // A = theta * (M * dxdot_dx + K) + dxdot_dx * _M_bar
        BLAS::copy(M, A);
        BLAS::aypx(A, dxdot_dx, K); // A = M * dxdot_dx + K

        BLAS::scale(A, theta); // A *= theta
        BLAS::axpy(A, dxdot_dx, _M_bar); // A += dxdot_dx * _M_bar
    }

    void getRhs(const Matrix &M, const Matrix &/*K*/, const Vector& b, Vector& rhs) const override
    {
        namespace BLAS = MathLib::BLAS;

        _crank_nicolson.getWeightedOldX(_tmp);

        auto const  theta          = _crank_nicolson.getTheta();

        // rhs = theta * (b + M * weighted_old_x) + _M_bar * weighted_old_x - _b_bar;
        BLAS::matMultAdd(M, _tmp, b, rhs); // rhs = b + M * weighted_old_x

        BLAS::scale(rhs, theta); // rhs *= theta
        BLAS::matMultAdd(_M_bar, _tmp, rhs, rhs); // rhs += _M_bar * weighted_old_x
        BLAS::axpy(rhs, -1.0, _b_bar); // rhs -= b
    }

    void getResidual(Matrix const& M, Matrix const& K, Vector const& b,
                     Vector const& x_new_timestep,
                     Vector& res) const override
    {
        namespace BLAS = MathLib::BLAS;

        auto const  alpha  = _crank_nicolson.getCurrentXWeight();
        auto const& x_curr = _crank_nicolson.getCurrentX(x_new_timestep);
        auto const  theta  = _crank_nicolson.getTheta();

        // x_dot  = alpha*x_new_timestep - x_old
        _crank_nicolson.getWeightedOldX(_tmp);
        BLAS::axpby(_tmp, alpha, -1.0, x_new_timestep);

        // res = theta * (M * x_dot + K*x_curr - b) + _M_bar * x_dot + _b_bar
        BLAS::matMult(M, _tmp, res); // res = M * x_dot; // the local vector x_dot seems to be necessary because of this multiplication
        BLAS::matMultAdd(K, x_curr, res, res); // res += K * x_curr
        BLAS::axpy(res, -1.0, b); // res = M * x_dot + K * x_curr - b

        BLAS::aypx(res, theta, _b_bar); // res = res * theta + _b_bar
        BLAS::matMultAdd(_M_bar, _tmp, res, res); // rs += _M_bar * x_dot
    }

    void getJacobian(Matrix const& Jac_in, Matrix& Jac_out) const override
    {
        namespace BLAS = MathLib::BLAS;

        auto const dxdot_dx = _crank_nicolson.getCurrentXWeight();
        auto const theta    = _crank_nicolson.getTheta();

        // J = theta * Jac + dxdot_dx * _M_bar
        BLAS::copy(Jac_in, Jac_out);
        BLAS::scale(Jac_out, theta);
        BLAS::axpy(Jac_out, dxdot_dx, _M_bar);
    }

    void pushMatrices(Matrix const& M, Matrix const& K, Vector const& b) override
    {
        namespace BLAS = MathLib::BLAS;

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
    CrankNicolson<Vector> const& _crank_nicolson;

    Matrix _M_bar;
    Vector _b_bar;
    mutable Vector _tmp; // used to store intermediate calculation results
};


template<typename Matrix, typename Vector, ODESystemTag ODETag>
std::unique_ptr<MatrixTranslator<Matrix, Vector, ODETag>>
createMatrixTranslator(TimeDiscretization<Vector> const& timeDisc)
{
    if (auto* fwd_euler = dynamic_cast<ForwardEuler<Vector> const*>(&timeDisc))
    {
        return std::unique_ptr<MatrixTranslator<Matrix, Vector, ODETag>>(
                new MatrixTranslatorForwardEuler<Matrix, Vector, ODETag>(*fwd_euler));
    }
    else if (auto* crank = dynamic_cast<CrankNicolson<Vector> const*>(&timeDisc))
    {
        return std::unique_ptr<MatrixTranslator<Matrix, Vector, ODETag>>(
                new MatrixTranslatorCrankNicolson<Matrix, Vector, ODETag>(*crank));
    }
    else
    {
        return std::unique_ptr<MatrixTranslator<Matrix, Vector, ODETag>>(
                new MatrixTranslatorGeneral<Matrix, Vector, ODETag>(timeDisc));
    }
}

//! @}

}
