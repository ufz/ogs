/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef NUMLIB_MATRIXTRANSLATOR_H
#define NUMLIB_MATRIXTRANSLATOR_H

#include <memory>

#include "Types.h"
#include "TimeDiscretization.h"

#include "MathLib/LinAlg/BLAS.h"

namespace NumLib
{

//! \addtogroup ODESolver
//! @{


/*! Translates matrices assembled by a provided ODE (or other equation)
 *  to some other matrices suitable to be passed on to nonlinear solvers.
 *
 * \tparam Matrix the type of matrices occuring in the linearization of the equation.
 * \tparam Vector the type of the solution vector of the ODEequation
 * \tparam ODETag a tag indicating the type of equation.
 */
template<typename Matrix, typename Vector, ODESystemTag ODETag>
class MatrixTranslator;


/*! Translates matrices assembled by a provided first order implicit quasi-linear ODE
 *  to some other matrices suitable to be passed on to nonlinear solvers.
 *
 * \tparam Matrix the type of matrices occuring in the linearization of the equation.
 * \tparam Vector the type of the solution vector of the ODEequation
 *
 * \see ODESystemTag::FirstOrderImplicitQuasilinear
 */
template<typename Matrix, typename Vector>
class MatrixTranslator<Matrix, Vector, ODESystemTag::FirstOrderImplicitQuasilinear>
{
public:
    //! Computes \c A from \c M and \c K.
    virtual void getA(Matrix const& M, Matrix const& K, Matrix& A) const = 0;

    //! Computes \c rhs from \c M, \c K and \c b.
    virtual void getRhs(const Matrix &M, const Matrix &K, const Vector& b, Vector& rhs) const = 0;

    /*! Computes \c res from \c M, \c K, \c b, \f$ \hat x \f$ and \f$ x_N \f$.
     * You might also want read the remarks on
     * \ref concept_time_discretization "time discretization".
     */
    virtual void getResidual(Matrix const& M, Matrix const& K, Vector const& b,
                             Vector const& x_new_timestep, Vector const& xdot,
                             Vector& res) const = 0;

    //! Computes the Jacobian of the residual and writes it to \c Jac_out.
    virtual void getJacobian(Matrix const& Jac_in, Matrix& Jac_out) const = 0;

    // TODO check (void)
    /*! Allows to store the given matrices internally for later use.
     *
     * \remark
     * This method has been provided in order to be able to implement the CrankNicolson
     * scheme.
     */
    virtual void pushMatrices(Matrix const& M, Matrix const& K, Vector const& b)
    {
        (void) M; (void) K; (void) b;
    }

    virtual ~MatrixTranslator() = default;
};


/*! General matrix translator used with time discretization schemes that have no special needs.
 *
 * \tparam Matrix the type of matrices occuring in the linearization of the equation.
 * \tparam Vector the type of the solution vector of the ODEequation
 * \tparam ODETag a tag indicating the type of equation.
 */
template<typename Matrix, typename Vector, ODESystemTag ODETag>
class MatrixTranslatorGeneral;


/*! General matrix translator for first order implicit quasi-linear ODEs, used with
 *  time discretization schemes that have no special needs.
 *
 * \tparam Matrix the type of matrices occuring in the linearization of the equation.
 * \tparam Vector the type of the solution vector of the ODEequation
 * \tparam ODETag a tag indicating the type of equation.
 *
 * \see ODESystemTag::FirstOrderImplicitQuasilinear
 * \remark
 * You might also want read the remarks on
 * \ref concept_time_discretization "time discretization".
 */
template<typename Matrix, typename Vector>
class MatrixTranslatorGeneral<Matrix, Vector, ODESystemTag::FirstOrderImplicitQuasilinear>
        : public MatrixTranslator<Matrix, Vector, ODESystemTag::FirstOrderImplicitQuasilinear>
{
public:
    /*! Constructs a new instance.
     *
     * \param timeDisc the time discretization scheme to be used.
     */
    MatrixTranslatorGeneral(TimeDiscretization<Vector> const& timeDisc)
        : _time_disc(timeDisc)
    {}

    //! Computes \f$ A = M \cdot \alpha + K \f$.
    void getA(Matrix const& M, Matrix const& K, Matrix& A) const override
    {
        namespace BLAS = MathLib::BLAS;

        auto const dxdot_dx = _time_disc.getCurrentXWeight();

        // A = M * dxdot_dx + K
        BLAS::copy(M, A);
        BLAS::aypx(A, dxdot_dx, K);
    }

    //! Computes \f$ \mathtt{rhs} = M \cdot x_O + b \f$.
    void getRhs(const Matrix &M, const Matrix &/*K*/, const Vector& b, Vector& rhs) const override
    {
        namespace BLAS = MathLib::BLAS;

        _time_disc.getWeightedOldX(_tmp);

        // rhs = M * weighted_old_x + b
        BLAS::matMultAdd(M, _tmp, b, rhs);
    }

    //! Computes \f$ r = M \cdot \hat x + K \cdot x_C - b \f$.
    void getResidual(Matrix const& M, Matrix const& K, Vector const& b,
                     Vector const& x_new_timestep,  Vector const& xdot,
                     Vector& res) const override
    {
        namespace BLAS = MathLib::BLAS;

        auto const& x_curr = _time_disc.getCurrentX(x_new_timestep);

        // res = M * x_dot + K * x_curr - b
        BLAS::matMult(M, xdot, res); // the local vector x_dot seems to be necessary because of this multiplication
        BLAS::matMultAdd(K, x_curr, res, res);
        BLAS::axpy(res, -1.0, b);
    }

    //! Writes \c Jac_in to \c Jac_out.
    //! \todo Do not copy.
    void getJacobian(Matrix const& Jac_in, Matrix& Jac_out) const override
    {
        namespace BLAS = MathLib::BLAS;

        BLAS::copy(Jac_in, Jac_out);
    }

private:
    TimeDiscretization<Vector> const& _time_disc; //!< the time discretization used.
    mutable Vector _tmp; //!< used to store intermediate calculation results
};


/*! Matrix translator used with the ForwardEuler scheme.
 *
 * \tparam Matrix the type of matrices occuring in the linearization of the equation.
 * \tparam Vector the type of the solution vector of the ODEequation
 * \tparam ODETag a tag indicating the type of equation.
 */
template<typename Matrix, typename Vector, ODESystemTag ODETag>
class MatrixTranslatorForwardEuler;


/*! Matrix translator for first order implicit quasi-linear ODEs,
 *  used with the ForwardEuler scheme.
 *
 * \tparam Matrix the type of matrices occuring in the linearization of the equation.
 * \tparam Vector the type of the solution vector of the ODEequation
 * \tparam ODETag a tag indicating the type of equation.
 *
 * \see ODESystemTag::FirstOrderImplicitQuasilinear
 * \remark
 * You might also want read the remarks on
 * \ref concept_time_discretization "time discretization".
 */
template<typename Matrix, typename Vector>
class MatrixTranslatorForwardEuler<Matrix, Vector, ODESystemTag::FirstOrderImplicitQuasilinear>
        : public MatrixTranslator<Matrix, Vector, ODESystemTag::FirstOrderImplicitQuasilinear>
{
public:
    /*! Constructs a new instance.
     *
     * \param timeDisc the time discretization scheme to be used.
     */
    MatrixTranslatorForwardEuler(ForwardEuler<Vector> const& timeDisc)
        : _fwd_euler(timeDisc)
    {}

    //! Computes \f$ A = M \cdot \alpha \f$.
    void getA(Matrix const& M, Matrix const& /*K*/, Matrix& A) const override
    {
        namespace BLAS = MathLib::BLAS;

        auto const dxdot_dx = _fwd_euler.getCurrentXWeight();

        // A = M * dxdot_dx
        BLAS::copy(M, A);
        BLAS::scale(A, dxdot_dx);
    }

    //! Computes \f$ \mathtt{rhs} = M \cdot x_O - K \cdot x_O + b \f$.
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

    //! Computes \f$ r = M \cdot \hat x + K \cdot x_C - b \f$.
    void getResidual(Matrix const& M, Matrix const& K, Vector const& b,
                     Vector const& x_new_timestep, Vector const& xdot,
                     Vector& res) const override
    {
        namespace BLAS = MathLib::BLAS;

        auto const& x_curr = _fwd_euler.getCurrentX(x_new_timestep);

        // res = M * x_dot + K * x_curr - b
        BLAS::matMult(M, xdot, res);
        BLAS::matMultAdd(K, x_curr, res, res);
        BLAS::axpy(res, -1.0, b);
    }

    //! Writes \c Jac_in to \c Jac_out.
    //! \todo Do not copy.
    void getJacobian(Matrix const& Jac_in, Matrix& Jac_out) const override
    {
        namespace BLAS = MathLib::BLAS;

        BLAS::copy(Jac_in, Jac_out);
    }

private:
    ForwardEuler<Vector> const& _fwd_euler; //!< the time discretization used.
    mutable Vector _tmp; //!< used to store intermediate calculation results
};


/*! Matrix translator used with the CrankNicolson scheme.
 *
 * \tparam Matrix the type of matrices occuring in the linearization of the equation.
 * \tparam Vector the type of the solution vector of the ODEequation
 * \tparam ODETag a tag indicating the type of equation.
 */
template<typename Matrix, typename Vector, ODESystemTag ODETag>
class MatrixTranslatorCrankNicolson;

/*! Matrix translator for first order implicit quasi-linear ODEs,
 *  used with the CrankNicolson scheme.
 *
 * \tparam Matrix the type of matrices occuring in the linearization of the equation.
 * \tparam Vector the type of the solution vector of the ODEequation
 * \tparam ODETag a tag indicating the type of equation.
 *
 * \see ODESystemTag::FirstOrderImplicitQuasilinear
 * \remark
 * You might also want read the remarks on
 * \ref concept_time_discretization "time discretization".
 */
template<typename Matrix, typename Vector>
class MatrixTranslatorCrankNicolson<Matrix, Vector, ODESystemTag::FirstOrderImplicitQuasilinear>
        : public MatrixTranslator<Matrix, Vector, ODESystemTag::FirstOrderImplicitQuasilinear>
{
public:
    /*! Constructs a new instance.
     *
     * \param timeDisc the time discretization scheme to be used.
     */
    MatrixTranslatorCrankNicolson(CrankNicolson<Vector> const& timeDisc)
        : _crank_nicolson(timeDisc)
    {}

    //! Computes \f$ A = \theta \cdot (M \cdot \alpha + K) + \bar M \cdot \alpha \f$.
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

    //! Computes \f$ \mathtt{rhs} = \theta \cdot (M \cdot x_O + b) + \bar M \cdot x_O - \bar b \f$.
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
    //! Computes \f$ r = \theta \cdot (M \cdot \hat x + K \cdot x_C - b) + \bar M \cdot \hat x + \bar b \f$.
    void getResidual(Matrix const& M, Matrix const& K, Vector const& b,
                     Vector const& x_new_timestep, Vector const& xdot,
                     Vector& res) const override
    {
        namespace BLAS = MathLib::BLAS;

        auto const& x_curr = _crank_nicolson.getCurrentX(x_new_timestep);
        auto const  theta  = _crank_nicolson.getTheta();

        // res = theta * (M * x_dot + K*x_curr - b) + _M_bar * x_dot + _b_bar
        BLAS::matMult(M, xdot, res); // res = M * x_dot
        BLAS::matMultAdd(K, x_curr, res, res); // res += K * x_curr
        BLAS::axpy(res, -1.0, b); // res = M * x_dot + K * x_curr - b

        BLAS::aypx(res, theta, _b_bar); // res = res * theta + _b_bar
        BLAS::matMultAdd(_M_bar, xdot, res, res); // rs += _M_bar * x_dot
    }

    /*! Computes \f$ \mathtt{Jac\_out} = \theta \cdot \mathtt{Jac\_in} + \bar M \cdot \alpha \f$.
     *
     * Where \c Jac_in is the Jacobian as assembled by the ODE system, i.e. in the same
     * fashion as for the BackwardEuler scheme.
     */
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

    /*! Saves internal state for use in the successive timestep;
     *  computes \f$ \bar M \f$ and \f$ \bar b \f$.
     *
     * \f$ \bar M \f$ and \f$ \bar b \f$ are computed as follows:
     *  \f{align}{
     *    \bar M &= (1-\theta) \cdot M \\
     *    \bar b &= (1-\theta) \cdot ( K \cdot x_n - b )
     *  \f}
     *
     * Where \f$ x_n \f$ is the solution at the timestep just finished.
     */
    void pushMatrices(Matrix const& M, Matrix const& K, Vector const& b) override
    {
        namespace BLAS = MathLib::BLAS;

        auto const theta = _crank_nicolson.getTheta();

        // Note: using x_old here is correct, since this method is called from within
        //       CrankNicolson::pushState() __after__ x_old has been updated to the result
        //       from the timestep just finished.
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

    Matrix _M_bar; //!< Used to adjust matrices and vectors assembled by the ODE.
                   //!< \see pushMatrices()
    Vector _b_bar; //!< Used to adjust vectors assembled by the ODE.
                   //!< \see pushMatrices()
    mutable Vector _tmp; //!< used to store intermediate calculation results
};


//! Creates a matrix translator suitable to work together with the given
//! time discretization scheme.
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

#endif // NUMLIB_MATRIXTRANSLATOR_H
