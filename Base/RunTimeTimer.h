#ifndef RUNTIMETIMER_H
#define RUNTIMETIMER_H

#include <sys/time.h>

#include "TimeMeasurementBase.h"

class RunTimeTimer
{
public:
        virtual void start();
        virtual void stop();
        virtual double elapsed();
	~RunTimeTimer() {};
private:
	timeval _start;
	timeval _stop;
};

#endif

