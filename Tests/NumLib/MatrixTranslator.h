#pragma once

#include <memory>

#include "ODETypes.h"
#include "TimeDiscretization.h"


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

    virtual Matrix getJacobian(Matrix Jac) const = 0;

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
