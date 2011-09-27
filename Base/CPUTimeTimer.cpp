#include "CPUTimeTimer.h"

void CPUTimeTimer::start()
{
	_start = clock();
}

void CPUTimeTimer::stop()
{
	_stop = clock();
}

double CPUTimeTimer::elapsed()
{
	return (_stop-_start)/(double)(CLOCKS_PER_SEC);
}

