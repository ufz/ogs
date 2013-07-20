/**
 * \file
 * \author Thomas Fischer
 * \date   2012-05-10
 * \brief  Definition of the RunTime class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef RUNTIME_H
#define RUNTIME_H

#ifndef _MSC_VER
#include <sys/time.h>
#else
#include <windows.h>
#endif

namespace BaseLib {

class RunTime
{
public:
	void start();
	void stop();
	double elapsed();
private:
#ifndef _MSC_VER
	timeval _start;
	timeval _stop;
#else
	unsigned long _start;
	unsigned long _stop;
#endif
};

} // end namespace BaseLib

#endif
