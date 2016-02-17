/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef NUMLIB_TYPES_H
#define NUMLIB_TYPES_H

namespace NumLib
{

//! \addtogroup ODESolver
//! @{

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
    FirstOrderImplicitQuasilinear,
    DUMMY, // TODO document
    NeumannBC // Sure, that's misuse of this enum, so sue me!
};

//! @}


// TODO move MatrixTraits to NumericsConfig.h?
template<typename Matrix>
struct MatrixTraits
/*
{
    using Index = int;
}
// */
;

} // namespace NumLib



#ifdef OGS_USE_EIGEN

#include<Eigen/Core>

namespace NumLib
{
template<>
struct MatrixTraits<Eigen::MatrixXd>
{
    using Index = Eigen::MatrixXd::Index;
};
}

#endif


#ifdef USE_PETSC

namespace MathLib { class PETScMatrix; }

namespace NumLib
{
template<>
struct MatrixTraits<MathLib::PETScMatrix>
{
    using Index = MathLib::PETScMatrix::IndexType;
};
}

#elif defined(OGS_USE_EIGEN)

namespace MathLib { class EigenMatrix; }

namespace NumLib
{
template<>
struct MatrixTraits<MathLib::EigenMatrix>
{
    using Index = MathLib::EigenMatrix::IndexType;
};
}

#endif


#endif // NUMLIB_TYPES_H
