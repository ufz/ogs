/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www./**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.net/LICENSE.txt
 *
 *
opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www./**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.net/LICENSE.txt
 *
 *
opengeosys.com/LICENSE.txt
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

