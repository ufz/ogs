/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef BASELIB_TMPUTIL_H
#define BASELIB_TMPUTIL_H

namespace BaseLib
{
//! Has sequence of integers as template parameters
template <int...>
struct IntegerSequence {
};

//! Generates an IntegerSequence.
//!
//! \see http://stackoverflow.com/a/7858971
template <int N, int... S>
struct GenerateIntegerSequence {
    // effectively pushes N-1 from the left to the list int... S of integers.
    typedef typename GenerateIntegerSequence<N - 1, N - 1, S...>::type type;
};

template <int... S>
struct GenerateIntegerSequence<0, S...> {
    typedef IntegerSequence<S...> type;
};
/* The template metaprogram proceeds in the following way:
 *
 * GenerateIntegerSequence<sizeof...(Args)>::type
 *
 * Assume sizeof...(Args) == 3. Let GIS := GenerateIntegerSequence
 * GIS<3, []>
 * -> GIS<2, [2]>
 * -> GIS<1, [1, 2]>
 * -> GIS<0, [0, 1, 2], which has member typedef IntegerSequence<0, 1, 2>
 */

}  // namespace BaseLib

#endif  // BASELIB_TMPUTIL_H
