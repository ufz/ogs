/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file RunTime.h
 *
 * Created on 2012-05-10 by Thomas Fischer
 */

#ifndef RUNTIME_H
#define RUNTIME_H

#include "TimeMeasurementBase.h"

#ifndef _WIN32
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
#ifndef _WIN32
	timeval _start;
	timeval _stop;
#else
	unsigned long _start;
	unsigned long _stop;
#endif
};

} // end namespace BaseLib

#endif
