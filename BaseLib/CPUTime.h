/**
 * \file
 * \author Thomas Fischer
 * \date   2012-05-10
 * \brief  Definition of the CPUTime class.
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef CPUTIME_H
#define CPUTIME_H

#include <ctime>

namespace BaseLib {

class CPUTime
{
public:
	void start();
    void stop();
    double elapsed();
private:
	clock_t _start;
	clock_t _stop;
};

} // end namespace BaseLib

#endif

