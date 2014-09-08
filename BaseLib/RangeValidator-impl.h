/**
 * \file
 * \author Thomas Fischer
 * \date   Sep 28, 2012
 * \brief  Implementation of the RangeValidator class.
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

namespace BaseLib {

template <typename NUMERIC_TYPE>
RangeValidator<NUMERIC_TYPE>::RangeValidator(NUMERIC_TYPE lower_limit, NUMERIC_TYPE upper_limit) :
	_lower_limit(lower_limit), _upper_limit(upper_limit)
{
}

template <typename NUMERIC_TYPE>
RangeValidator<NUMERIC_TYPE>::~RangeValidator()
{}

template <typename NUMERIC_TYPE>
void RangeValidator<NUMERIC_TYPE>::resetLowerLimits(NUMERIC_TYPE lower_limit)
{
	_lower_limit = lower_limit;
}

template <typename NUMERIC_TYPE>
void RangeValidator<NUMERIC_TYPE>:: resetUpperLimits(NUMERIC_TYPE upper_limit)
{
	_upper_limit = upper_limit;
}

template <typename NUMERIC_TYPE>
bool RangeValidator<NUMERIC_TYPE>::isValidValue (NUMERIC_TYPE test_value) const
{
	if (test_value < _lower_limit)
		return false;

	if (_upper_limit < test_value)
		return false;

	return true;
}

} // end namespace BaseLib
