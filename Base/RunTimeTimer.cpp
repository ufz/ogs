#include "RunTimeTimer.h"

void RunTimeTimer::start()
{
#ifndef _WIN32
	gettimeofday(&_start, 0);
#else
	_start = timeGetTime();
#endif
}

void RunTimeTimer::stop()
{
#ifndef _WIN32
	gettimeofday(&_stop, 0);
#else
	_stop = timeGetTime();
#endif
}

double RunTimeTimer::elapsed()
{
#ifndef _WIN32
	return (_stop.tv_sec + _stop.tv_usec/1000000.0 - (_start.tv_sec + _start.tv_usec/1000000.0));
#else
	return (_stop - _start) / 1000;
#endif
}
