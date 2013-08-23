/**
 * \author Norihiro Watanabe
 * \date   2013-04-16
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef ASSEMBLERLIB_SERIALEXECUTOR_H_H
#define ASSEMBLERLIB_SERIALEXECUTOR_H_H

#include <vector>

namespace AssemblerLib
{

/// Executes a \c f for each element in input vector.
/// Return values of the function call are ignored.
///
/// \tparam T type of input vector elements
/// \tparam F \c f type
///
/// \param v a vector of T pointers
/// \param f a function that accepts a pointer to T and an index as arguments
template <typename T, typename F>
void
serialExecute(std::vector<T*> const& v, F const& f)
{
    for (std::size_t i = 0; i < v.size(); i++)
        f(v[i], i);
};

}   // namespace AssemblerLib

#endif  // ASSEMBLERLIB_SERIALEXECUTOR_H_H
