/**
 * \author Norihiro Watanabe
 * \date   2012-08-03
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#ifndef ITIMESTEPALGORITHM_H_
#define ITIMESTEPALGORITHM_H_

#include <vector>

#include "NumLib/TimeStepping/TimeStep.h"

namespace NumLib
{

/**
 * \brief Interface of time stepping algorithms
 *
 */
class ITimeStepAlgorithm
{
public:
    /// return the beginning of time steps
    virtual double begin() const = 0;

    /// return the end of time steps
    virtual double end() const = 0;

    /// move to the next time step
    /// \return true if the next step exists
    virtual bool next() = 0;

    /// return current time step
    virtual const TimeStep getTimeStep() const = 0;

    /// return if current time step is accepted or not
    virtual bool accepted() const = 0;

    /// return a history of time step sizes
    virtual const std::vector<double>& getTimeStepSizeHistory() const = 0;

    virtual ~ITimeStepAlgorithm() {}
};

} //NumLib

#endif //ITIMESTEPALGORITHM_H_
