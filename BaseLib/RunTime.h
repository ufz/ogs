/**
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 * \file RunTime.h
 *
 * Created on 2012-05-10 by Thomas Fischer
 */

#ifndef RUNTIME_H
#define RUNTIME_H

#include "TimeMeasurementBase.h"

#ifndef _MSC_VER
#include <sys/time.h>
#else
#include <windows.h>
#endif

#include "TimeMeasurementBase.h"

namespace BaseLib {

class RunTime : public TimeMeasurementBase
{
public:
	virtual void start();
	virtual void stop();
	virtual double elapsed();
	~RunTime() {};
private:
#ifndef _MSC_VER
	timeval _start;
	timeval _stop;
#else
	unsigned long _start;
	unsigned long _stop;
#endif
};

} // end namespace BaseLib

#endif
