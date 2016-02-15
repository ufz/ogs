#pragma once

namespace NumLib
{

//! \addtogroup ODESolver
//! @{

// TODO subject to change
using IndexType = std::size_t;

//! Tag used to specify which nonlinear solver will be used.
enum class NonlinearSolverTag : bool {
    Picard /*!< Picard fixpoint iteration scheme */,
    Newton /*!< Newton-Raphson iteration scheme */
};

//! Tag used to specify the type of ODE.
enum class ODESystemTag : char
{
    /*! First order implicit quasi-linear ODE
     *
     * This is an ODE of the form
     * \f$ M(x,t)\cdot \dot x + K(x,t) \cdot x - b(x,t)
     *  =: r(\dot x, x, t) \stackrel{!}{=} 0 \f$
     */
    FirstOrderImplicitQuasilinear
};

//! @}

}
