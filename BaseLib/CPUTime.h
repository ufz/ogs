/**
 * \file
 * \author Thomas Fischer
 * \date   2012-05-10
 * \brief  Definition of the CPUTime class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef CPUTIME_H
#define CPUTIME_H

#include <ctime>

#include "TimeMeasurementBase.h"

namespace BaseLib {

class CPUTime
{
public:
	virtual void start();
    virtual void stop();
    virtual double elapsed();
	~CPUTime() {};
private:
	clock_t _start;
	clock_t _stop;
};

} // end namespace BaseLib

#endif

