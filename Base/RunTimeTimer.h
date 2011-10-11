#ifndef RUNTIMETIMER_H
#define RUNTIMETIMER_H

#include "TimeMeasurementBase.h"

#ifndef _WIN32
#include <sys/time.h>
#else
#include <windows.h>
#endif

#include "TimeMeasurementBase.h"

class RunTimeTimer : public TimeMeasurementBase
{
public:
	virtual void start();
	virtual void stop();
	virtual double elapsed();
	~RunTimeTimer() {};
private:
#ifndef _WIN32
	timeval _start;
	timeval _stop;
#else
	unsigned long _start;
	unsigned long _stop;
#endif
};

#endif
