/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 * \file CPUTime.h
 *
 * Created on 2012-05-10 by Thomas Fischer
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

