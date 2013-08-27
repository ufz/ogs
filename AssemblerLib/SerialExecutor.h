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

namespace AssemblerLib
{

struct SerialExecutor
{
    /// Executes a \c f for each element from the input container.
    /// Return values of the function call are ignored.
    ///
    /// \tparam F   \c f type.
    /// \tparam C   input container type.
    ///
    /// \param f    a function that accepts a pointer to container's elements and
    ///             an index as arguments.
    /// \param c    a container supporting access over operator[].
    template <typename F, typename C>
    static
    void
    execute(F const& f, C const& c)
    {
        for (std::size_t i = 0; i < c.size(); i++)
            f(c[i], i);
    };

};

}   // namespace AssemblerLib

#endif  // ASSEMBLERLIB_SERIALEXECUTOR_H_H
