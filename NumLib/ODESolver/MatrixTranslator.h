/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <memory>

#include "TimeDiscretization.h"
#include "Types.h"

#include "MathLib/LinAlg/LinAlg.h"

namespace NumLib
{
//! \addtogroup ODESolver
//! @{

/*! Translates matrices and vectors assembled by a provided ODE (or other
 * equation) to the stiffness matrix and right-hand side vector of a linear
 * equation system that can be solved by a linear equation system solver.
 *
 * \tparam ODETag a tag indicating the type of equation.
 */
template <ODESystemTag ODETag>
class MatrixTranslator;

/*! Translates matrices assembled by a provided first order implicit
 * quasi-linear ODE to some other matrices suitable to be passed on to nonlinear solvers.
 *
 * \see ODESystemTag::FirstOrderImplicitQuasilinear
 */
template <>
class MatrixTranslator<ODESystemTag::FirstOrderImplicitQuasilinear>
{
public:
    //! Computes \c A from \c M and \c K.
    virtual void computeA(GlobalMatrix const& M, GlobalMatrix const& K,
                          GlobalMatrix& A) const = 0;

    //! Computes \c rhs from \c M, \c K and \c b.
    virtual void computeRhs(const GlobalMatrix& M, const GlobalMatrix& K,
                            const GlobalVector& b, GlobalVector& rhs) const = 0;

    /*! Computes \c res from \c M, \c K, \c b, \f$ \hat x \f$ and \f$ x_N \f$.
     * You might also want read the remarks on
     * \ref concept_time_discretization "time discretization".
     */
    virtual void computeResidual(GlobalMatrix const& M, GlobalMatrix const& K,
                                 GlobalVector const& b,
                                 GlobalVector const& x_new_timestep,
                                 GlobalVector const& xdot,
                                 GlobalVector& res) const = 0;

    //! Computes the Jacobian of the residual and writes it to \c Jac_out.
    virtual void computeJacobian(GlobalMatrix const& Jac_in,
                                 GlobalMatrix& Jac_out) const = 0;

    /*! Allows to store the given matrices internally for later use.
     *
     * \remark
     * This method has been provided in order to be able to implement the
     * CrankNicolson
     * scheme.
     */
    virtual void pushMatrices(GlobalMatrix const& /*M*/,
                              GlobalMatrix const& /*K*/,
                              GlobalVector const& /*b*/)
    {
    }

    virtual ~MatrixTranslator() = default;
};

/*! General GlobalMatrix translator used with time discretization schemes that
 * have no special needs.
 *
 * \tparam ODETag a tag indicating the type of equation.
 */
template <ODESystemTag ODETag>
class MatrixTranslatorGeneral;

/*! General matrix translator for first order implicit quasi-linear ODEs, used
 * with time discretization schemes that have no special needs.
 *
 * \see ODESystemTag::FirstOrderImplicitQuasilinear
 * \remark
 * You might also want read the remarks on
 * \ref concept_time_discretization "time discretization".
 */
template <>
class MatrixTranslatorGeneral<ODESystemTag::FirstOrderImplicitQuasilinear>
    : public MatrixTranslator<ODESystemTag::FirstOrderImplicitQuasilinear>
{
public:
    /*! Constructs a new instance.
     *
     * \param timeDisc the time discretization scheme to be used.
     */
    MatrixTranslatorGeneral(TimeDiscretization const& timeDisc)
        : _time_disc(timeDisc)
    {
    }

    //! Computes \f$ A = M \cdot \alpha + K \f$.
    void computeA(GlobalMatrix const& M, GlobalMatrix const& K,
                  GlobalMatrix& A) const override;

    //! Computes \f$ \mathtt{rhs} = M \cdot x_O + b \f$.
    void computeRhs(const GlobalMatrix& M, const GlobalMatrix& /*K*/,
                    const GlobalVector& b, GlobalVector& rhs) const override;

    //! Computes \f$ r = M \cdot \hat x + K \cdot x_C - b \f$.
    void computeResidual(GlobalMatrix const& M, GlobalMatrix const& K,
                         GlobalVector const& b,
                         GlobalVector const& x_new_timestep,
                         GlobalVector const& xdot,
                         GlobalVector& res) const override;

    //! Writes \c Jac_in to \c Jac_out.
    //! \todo Do not copy.
    void computeJacobian(GlobalMatrix const& Jac_in,
                         GlobalMatrix& Jac_out) const override;

private:
    TimeDiscretization const&
        _time_disc;  //!< the time discretization used.

    //! ID of the vector storing intermediate computations.
    mutable std::size_t _tmp_id = 0u;
};

/*! GlobalMatrix translator used with the ForwardEuler scheme.
 *
 * \tparam ODETag a tag indicating the type of equation.
 */
template <ODESystemTag ODETag>
class MatrixTranslatorForwardEuler;

/*! GlobalMatrix translator for first order implicit quasi-linear ODEs,
 *  used with the ForwardEuler scheme.
 *
 * \see ODESystemTag::FirstOrderImplicitQuasilinear
 * \remark
 * You might also want read the remarks on
 * \ref concept_time_discretization "time discretization".
 */
template <>
class MatrixTranslatorForwardEuler<ODESystemTag::FirstOrderImplicitQuasilinear>
    : public MatrixTranslator<ODESystemTag::FirstOrderImplicitQuasilinear>
{
public:
    /*! Constructs a new instance.
     *
     * \param timeDisc the time discretization scheme to be used.
     */
    MatrixTranslatorForwardEuler(ForwardEuler const& timeDisc)
        : _fwd_euler(timeDisc)
    {
    }

    //! Computes \f$ A = M \cdot \alpha \f$.
    void computeA(GlobalMatrix const& M, GlobalMatrix const& /*K*/,
                  GlobalMatrix& A) const override;

    //! Computes \f$ \mathtt{rhs} = M \cdot x_O - K \cdot x_O + b \f$.
    void computeRhs(const GlobalMatrix& M, const GlobalMatrix& K,
                    const GlobalVector& b, GlobalVector& rhs) const override;

    //! Computes \f$ r = M \cdot \hat x + K \cdot x_C - b \f$.
    void computeResidual(GlobalMatrix const& M, GlobalMatrix const& K,
                         GlobalVector const& b,
                         GlobalVector const& x_new_timestep,
                         GlobalVector const& xdot,
                         GlobalVector& res) const override;

    //! Writes \c Jac_in to \c Jac_out.
    //! \todo Do not copy.
    void computeJacobian(GlobalMatrix const& Jac_in,
                         GlobalMatrix& Jac_out) const override;

private:
    ForwardEuler const& _fwd_euler;  //!< the time discretization used.

    //! ID of the vector storing intermediate computations.
    mutable std::size_t _tmp_id = 0u;
};

/*! GlobalMatrix translator used with the CrankNicolson scheme.
 *
 * \tparam ODETag a tag indicating the type of equation.
 */
template <ODESystemTag ODETag>
class MatrixTranslatorCrankNicolson;

/*! GlobalMatrix translator for first order implicit quasi-linear ODEs,
 *  used with the CrankNicolson scheme.
 *
 * \see ODESystemTag::FirstOrderImplicitQuasilinear
 * \remark
 * You might also want read the remarks on
 * \ref concept_time_discretization "time discretization".
 */
template <>
class MatrixTranslatorCrankNicolson<ODESystemTag::FirstOrderImplicitQuasilinear>
    : public MatrixTranslator<ODESystemTag::FirstOrderImplicitQuasilinear>
{
public:
    /*! Constructs a new instance.
     *
     * \param timeDisc the time discretization scheme to be used.
     */
    MatrixTranslatorCrankNicolson(CrankNicolson const& timeDisc)
        : _crank_nicolson(timeDisc),
          _M_bar(NumLib::GlobalMatrixProvider::provider
                     .getMatrix()),
          _b_bar(
              NumLib::GlobalVectorProvider::provider.getVector())
    {
    }

    ~MatrixTranslatorCrankNicolson() override
    {
        NumLib::GlobalMatrixProvider::provider.releaseMatrix(
            _M_bar);
        NumLib::GlobalVectorProvider::provider.releaseVector(
            _b_bar);
    }

    //! Computes \f$ A = \theta \cdot (M \cdot \alpha + K) + \bar M \cdot \alpha
    //! \f$.
    void computeA(GlobalMatrix const& M, GlobalMatrix const& K,
                  GlobalMatrix& A) const override;

    //! Computes \f$ \mathtt{rhs} = \theta \cdot (M \cdot x_O + b) + \bar M
    //! \cdot x_O - \bar b \f$.
    void computeRhs(const GlobalMatrix& M, const GlobalMatrix& /*K*/,
                    const GlobalVector& b, GlobalVector& rhs) const override;

    //! Computes \f$ r = \theta \cdot (M \cdot \hat x + K \cdot x_C - b) + \bar
    //! M \cdot \hat x + \bar b \f$.
    void computeResidual(GlobalMatrix const& M, GlobalMatrix const& K,
                         GlobalVector const& b,
                         GlobalVector const& x_new_timestep,
                         GlobalVector const& xdot,
                         GlobalVector& res) const override;

    /*! Computes \f$ \mathtt{Jac\_out} = \theta \cdot \mathtt{Jac\_in} + \bar M
     * \cdot \alpha \f$.
     *
     * Where \c Jac_in is the Jacobian as assembled by the ODE system, i.e. in
     * the same fashion as for the BackwardEuler scheme.
     */
    void computeJacobian(GlobalMatrix const& Jac_in,
                         GlobalMatrix& Jac_out) const override;

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
    void pushMatrices(GlobalMatrix const& M, GlobalMatrix const& K,
                      GlobalVector const& b) override;

private:
    CrankNicolson const& _crank_nicolson;

    GlobalMatrix&
        _M_bar;  //!< Used to adjust matrices and vectors assembled by the ODE.
                 //!< \see pushMatrices()
    GlobalVector& _b_bar;  //!< Used to adjust vectors assembled by the ODE.
                           //!< \see pushMatrices()

    //! ID of the vector storing intermediate computations.
    mutable std::size_t _tmp_id = 0u;
};

//! Creates a GlobalMatrix translator suitable to work together with the given
//! time discretization scheme.
template <ODESystemTag ODETag>
std::unique_ptr<MatrixTranslator<ODETag>> createMatrixTranslator(
    TimeDiscretization const& timeDisc)
{
    if (auto* fwd_euler = dynamic_cast<ForwardEuler const*>(&timeDisc))
    {
        return std::unique_ptr<MatrixTranslator<ODETag>>(
            new MatrixTranslatorForwardEuler<ODETag>(*fwd_euler));
    }
    else if (auto* crank =
                 dynamic_cast<CrankNicolson const*>(&timeDisc))
    {
        return std::unique_ptr<MatrixTranslator<ODETag>>(
            new MatrixTranslatorCrankNicolson<ODETag>(*crank));
    }
    else
    {
        return std::unique_ptr<MatrixTranslator<ODETag>>(
            new MatrixTranslatorGeneral<ODETag>(timeDisc));
    }
}

//! @}
}
