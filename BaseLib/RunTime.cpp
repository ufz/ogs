/**
 * \file
 * \author Thomas Fischer
 * \date   2012-05-10
 * \brief  Implementation of the RunTime class.
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "RunTime.h"

namespace BaseLib {

void RunTime::start()
{
#ifndef _MSC_VER
	gettimeofday(&_start, 0);
#else
	_start = timeGetTime();
#endif
}

void RunTime::stop()
{
#ifndef _MSC_VER
	gettimeofday(&_stop, 0);
#else
	_stop = timeGetTime();
#endif
}

double RunTime::elapsed()
{
#ifndef _MSC_VER
	return (_stop.tv_sec + _stop.tv_usec/1000000.0 - (_start.tv_sec + _start.tv_usec/1000000.0));
#else
	return (_stop - _start) / 1000.0;
#endif
}

} // end namespace BaseLib
