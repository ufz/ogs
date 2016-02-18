/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef NUMLIB_EQUATIONSYSTEM_H
#define NUMLIB_EQUATIONSYSTEM_H

namespace NumLib
{

//! \addtogroup ODESolver
//! @{

//! Collection of basic methods every equation system must provide.
class EquationSystem
{
public:
    //! Return the number of equations.
    virtual std::size_t getNumEquations() const = 0;

    /*! Check whether this is actually a linear equation system.
     *
     * \remark
     * Depending on its parameters an in general nonlinear equation system
     * can be linear in special cases. With this method it is possible to
     * detect that at runtime and thus save some computations.
     */
    virtual bool isLinear() const = 0;

    virtual ~EquationSystem() = default;
};

//! @}

}

#endif
