/**
 * \file
 * \author Thomas Fischer
 * \date   no date
 * \brief  Definition of the TimeMeasurementBase class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
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
