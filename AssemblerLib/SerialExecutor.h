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

/// Executes a \c f for each element from the input container.
/// Return values of the function call are ignored.
///
/// \tparam C   input container type.
/// \tparam F   \c f type.
///
/// \param c    a container supporting access over operator[].
/// \param f    a function that accepts a pointer to container's elements and
///             an index as arguments.
template <typename C, typename F>
void
serialExecute(C const& c, F const& f)
{
    for (std::size_t i = 0; i < c.size(); i++)
        f(c[i], i);
};

/// Executes a \c f for each element from the input range [first, last)
/// Return values of the function call are ignored.
///
/// \tparam I   input iterator type.
/// \tparam F   \c f type.
///
/// \param first    iterator to the first element.
/// \param last     iterator to one after the last element.
/// \param f        a function accepting a pointer to container's elements and
///                 an index as arguments.
template <typename I, typename F>
void
serialExecute(I first, I last, F const& f)
{
    std::size_t count = 0;
    while (first != last)
        f(*first++, count++);
};
}   // namespace AssemblerLib

#endif  // ASSEMBLERLIB_SERIALEXECUTOR_H_H
