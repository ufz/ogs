/*
 * CPUTime.h
 *
 *  Created on: May 10, 2012
 *      Author: TF
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

