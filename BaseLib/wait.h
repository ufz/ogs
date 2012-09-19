/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file wait.h
 *
 * Created on 2011-02-17 by Karsten Rink
 */


#ifndef WAIT_H
#define WAIT_H

#include <ctime>

namespace BaseLib {

void wait(int seconds)
{
	time_t start_time, cur_time;

	time(&start_time);
	do
	{
		 time(&cur_time);
	}
	while((cur_time - start_time) < seconds);
}

} // end namespace BaseLib

#endif //WAIT_H
