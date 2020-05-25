/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <memory>

#include "MatrixTranslator.h"
#include "NonlinearSystem.h"
#include "ODESystem.h"
#include "TimeDiscretization.h"

namespace NumLib
{
//! \addtogroup ODESolver
//! @{

/*! A NonlinearSystem together with some TimeDiscretization scheme.
 *
 * This is the interface of an ODE towards the TimeLoop.
 * This interface is abstract, it represents any type of first-order ODE.
 *
 * \tparam NLTag  a tag indicating the method used for resolving nonlinearities.
 */
template <NonlinearSolverTag NLTag>
class TimeDiscretizedODESystemBase : public NonlinearSystem<NLTag>
{
public:
    //! Exposes the used time discretization scheme.
    virtual TimeDiscretization& getTimeDiscretization() = 0;
};

/*! A NonlinearSystem together with some TimeDiscretization scheme.
 *
 * This class represents a specific type of first-order ODE, as indicated by
 * \c ODETag.
 *
 * \tparam ODETag a tag indicating the type of ODE.
 * \tparam NLTag  a tag indicating the method used for resolving nonlinearities.
 */
template <ODESystemTag ODETag, NonlinearSolverTag NLTag>
class TimeDiscretizedODESystem;

/*! Time discretized first order implicit quasi-linear ODE;
 *  to be solved using the Newton-Raphson method for resolving nonlinearities.
 *
 * \see ODESystemTag::FirstOrderImplicitQuasilinear
 */
template <>
class TimeDiscretizedODESystem<ODESystemTag::FirstOrderImplicitQuasilinear,
                               NonlinearSolverTag::Newton>
    final : public TimeDiscretizedODESystemBase<NonlinearSolverTag::Newton>
{
public:
    //! A tag indicating the type of ODE.
    static const ODESystemTag ODETag =
        ODESystemTag::FirstOrderImplicitQuasilinear;

    //! The type of ODE.
    using ODE = ODESystem<ODETag, NonlinearSolverTag::Newton>;
    //! The auxiliary class that computes the matrix/vector used by the
    //! nonlinear solver.
    using MatTrans = MatrixTranslator<ODETag>;
    //! A shortcut for a general time discretization scheme
    using TimeDisc = TimeDiscretization;

    /*! Constructs a new instance.
     *
     * \param process_id ID of the ODE to be solved.
     * \param ode the ODE to be wrapped.
     * \param time_discretization the time discretization to be used.
     */
    explicit TimeDiscretizedODESystem(const int process_id, ODE& ode,
                                      TimeDisc& time_discretization);

    ~TimeDiscretizedODESystem() override;

    void assemble(std::vector<GlobalVector*> const& x_new_timestep,
                  std::vector<GlobalVector*> const& x_prev,
                  int const process_id) override;

    void getResidual(GlobalVector const& x_new_timestep,
                     GlobalVector const& x_prev,
                     GlobalVector& res) const override;

    void getJacobian(GlobalMatrix& Jac) const override;

    void computeKnownSolutions(GlobalVector const& x,
                               int const process_id) override;

    void applyKnownSolutions(GlobalVector& x) const override;

    void applyKnownSolutionsNewton(GlobalMatrix& Jac, GlobalVector& res,
                                   GlobalVector& minus_delta_x) const override;

    bool isLinear() const override { return ode_.isLinear(); }

    void preIteration(const unsigned iter, GlobalVector const& x) override
    {
        ode_.preIteration(iter, x);
    }

    IterationResult postIteration(GlobalVector const& x) override
    {
        return ode_.postIteration(x);
    }

    TimeDisc& getTimeDiscretization() override { return time_disc_; }
    MathLib::MatrixSpecifications getMatrixSpecifications(
        const int process_id) const override
    {
        return ode_.getMatrixSpecifications(process_id);
    }

private:
    ODE& ode_;             //!< ode the ODE being wrapped
    TimeDisc& time_disc_;  //!< the time discretization to being used

    //! the object used to compute the matrix/vector for the nonlinear solver
    std::unique_ptr<MatTrans> mat_trans_;

    using Index = MathLib::MatrixVectorTraits<GlobalMatrix>::Index;
    std::vector<NumLib::IndexValueVector<Index>> const* known_solutions_ =
        nullptr;  //!< stores precomputed values for known solutions

    GlobalMatrix* Jac_;  //!< the Jacobian of the residual
    GlobalMatrix* M_;    //!< Matrix \f$ M \f$.
    GlobalMatrix* K_;    //!< Matrix \f$ K \f$.
    GlobalVector* b_;    //!< Matrix \f$ b \f$.

    std::size_t Jac_id_ = 0u;  //!< ID of the \c Jac_ matrix.
    std::size_t M_id_ = 0u;    //!< ID of the \c M_ matrix.
    std::size_t K_id_ = 0u;    //!< ID of the \c K_ matrix.
    std::size_t b_id_ = 0u;    //!< ID of the \c b_ vector.

    //! ID of the vector storing xdot in intermediate computations.
    mutable std::size_t xdot_id_ = 0u;
};

/*! Time discretized first order implicit quasi-linear ODE;
 *  to be solved using the Picard fixpoint iteration method for resolving
 * nonlinearities.
 *
 * \see ODESystemTag::FirstOrderImplicitQuasilinear
 */
template <>
class TimeDiscretizedODESystem<ODESystemTag::FirstOrderImplicitQuasilinear,
                               NonlinearSolverTag::Picard>
    final : public TimeDiscretizedODESystemBase<NonlinearSolverTag::Picard>
{
public:
    //! A tag indicating the type of ODE.
    static const ODESystemTag ODETag =
        ODESystemTag::FirstOrderImplicitQuasilinear;

    //! The type of ODE.
    using ODE = ODESystem<ODETag, NonlinearSolverTag::Picard>;
    //! The auxiliary class that computes the matrix/vector used by the
    //! nonlinear solver.
    using MatTrans = MatrixTranslator<ODETag>;
    //! A shortcut for a general time discretization scheme
    using TimeDisc = TimeDiscretization;

    //! Constructs a new instance.
    explicit TimeDiscretizedODESystem(const int process_id, ODE& ode,
                                      TimeDisc& time_discretization);

    ~TimeDiscretizedODESystem() override;

    void assemble(std::vector<GlobalVector*> const& x_new_timestep,
                  std::vector<GlobalVector*> const& x_prev,
                  int const process_id) override;

    void getA(GlobalMatrix& A) const override
    {
        mat_trans_->computeA(*M_, *K_, A);
    }

    void getRhs(GlobalVector const& x_prev, GlobalVector& rhs) const override
    {
        mat_trans_->computeRhs(*M_, *K_, *b_, x_prev, rhs);
    }

    void computeKnownSolutions(GlobalVector const& x,
                               int const process_id) override;

    void applyKnownSolutions(GlobalVector& x) const override;

    void applyKnownSolutionsPicard(GlobalMatrix& A, GlobalVector& rhs,
                                   GlobalVector& x) const override;

    bool isLinear() const override { return ode_.isLinear(); }

    void preIteration(const unsigned iter, GlobalVector const& x) override
    {
        ode_.preIteration(iter, x);
    }

    IterationResult postIteration(GlobalVector const& x) override
    {
        return ode_.postIteration(x);
    }

    TimeDisc& getTimeDiscretization() override { return time_disc_; }
    MathLib::MatrixSpecifications getMatrixSpecifications(
        const int process_id) const override
    {
        return ode_.getMatrixSpecifications(process_id);
    }

private:
    ODE& ode_;             //!< ode the ODE being wrapped
    TimeDisc& time_disc_;  //!< the time discretization to being used

    //! the object used to compute the matrix/vector for the nonlinear solver
    std::unique_ptr<MatTrans> mat_trans_;

    using Index = MathLib::MatrixVectorTraits<GlobalMatrix>::Index;
    std::vector<NumLib::IndexValueVector<Index>> const* known_solutions_ =
        nullptr;  //!< stores precomputed values for known solutions

    GlobalMatrix* M_;  //!< Matrix \f$ M \f$.
    GlobalMatrix* K_;  //!< Matrix \f$ K \f$.
    GlobalVector* b_;  //!< Matrix \f$ b \f$.

    std::size_t M_id_ = 0u;  //!< ID of the \c M_ matrix.
    std::size_t K_id_ = 0u;  //!< ID of the \c K_ matrix.
    std::size_t b_id_ = 0u;  //!< ID of the \c b_ vector.
};

//! @}
}  // namespace NumLib
