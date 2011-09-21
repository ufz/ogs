#include "RunTimeTimer.h"

void RunTimeTimer::start()
{
	gettimeofday(&_start, 0);
}

void RunTimeTimer::stop()
{
	gettimeofday(&_stop, 0);
}

double RunTimeTimer::elapsed()
{
	return (_stop.tv_sec + _stop.tv_usec/1000000.0 - (_start.tv_sec + _start.tv_usec/1000000.0));
}

