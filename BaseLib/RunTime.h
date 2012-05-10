/*
 * RunTime.h
 *
 *  Created on: May 10, 2012
 *      Author: TF
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
