/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file TimeMeasurementBase.h
 *
 * Created on xxxx-xx-xx by Thomas Fischer
 */

#ifndef TIMEMEASUREMENT_H
#define TIMEMEASUREMENT_H

namespace BaseLib {

class TimeMeasurementBase
{
public:
	virtual void start () = 0;
	virtual void stop () = 0;
	virtual double elapsed () = 0;
	virtual ~TimeMeasurementBase () {};
};

} // end namespace BaseLib

#endif
