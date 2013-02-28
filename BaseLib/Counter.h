/**
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef BASELIB_GLOBALCOUNTER_H
#define BASELIB_GLOBALCOUNTER_H

#include <cstddef>

namespace BaseLib
{
template <typename X>
struct Counter
{
	Counter()
	{
		_counter_value++;
	}
	static std::size_t _counter_value;
};

template <typename X> std::size_t Counter<X>::_counter_value(0);

} // end namespace BaseLib

#endif  // BASELIB_GLOBALCOUNTER_H
