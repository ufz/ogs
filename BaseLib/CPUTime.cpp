/**
 * \file
 * \author Thomas Fischer
 * \date   2012-05-10
 * \brief  Implementation of the CPUTime class.
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
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
