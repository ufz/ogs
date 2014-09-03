/**
 * \file
 * \author Thomas Fischer
 * \date   Sep 28, 2012
 * \brief  Definition of the RangeValidator calss.
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef RANGEVALIDATOR_H_
#define RANGEVALIDATOR_H_

namespace BaseLib {

template <typename NUMERIC_TYPE>
class RangeValidator {
public:
	RangeValidator(NUMERIC_TYPE lower_limit, NUMERIC_TYPE upper_limit);
	void resetLowerLimits(NUMERIC_TYPE lower_limit);
	void resetUpperLimits(NUMERIC_TYPE upper_limit);
	bool isValidValue (NUMERIC_TYPE test_value) const;
	virtual ~RangeValidator();
private:
	NUMERIC_TYPE _lower_limit;
	NUMERIC_TYPE _upper_limit;
};

}

#include "RangeValidator-impl.h"

#endif /* RANGEVALIDATOR_H_ */
