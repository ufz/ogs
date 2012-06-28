/**
 * \file CPUTime.cpp
 *
 * Created on 2012-05-10 by Thomas Fischer
 */

#include "CPUTime.h"

namespace BaseLib {

void CPUTime::start()
{
	_start = clock();
}

void CPUTime::stop()
{
	_stop = clock();
}

double CPUTime::elapsed()
{
	return (_stop-_start)/(double)(CLOCKS_PER_SEC);
}

} // end namespace BaseLib
