/**
 * \author Norihiro Watanabe
 * \date   2013-04-16
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
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
    template <typename F, typename C, typename ...Args_>
    static
    void
#if defined(_MSC_VER) && (_MSC_VER >= 1700)
    execute(F& f, C const& c, Args_&&... args)
#else
    execute(F const& f, C const& c, Args_&&... args)
#endif
    {
        for (std::size_t i = 0; i < c.size(); i++)
            f(i, c[i], std::forward<Args_>(args)...);
    }

    /// Same as execute(f, c), but with two containers, where the second one is
    /// modified.
    ///
    /// \tparam F    \c f type.
    /// \tparam C    input container type.
    /// \tparam Data input/output container type.
    ///
    /// \param f     a function that accepts a pointer to container's elements,
    ///              an index, and a second container element as arguments, which
    ///              is modified.
    /// \param c     a container supporting const access over operator[] and size().
    /// \param data  a container supporting non-const access over operator[] and size().
    template <typename F, typename C, typename Data, typename ...Args_>
    static
    void
#if defined(_MSC_VER) && (_MSC_VER >= 1700)
    execute(F& f, C const& c, Data& data, Args_&&... args)
#else
    execute(F const& f, C const& c, Data& data, Args_&&... args)
#endif
    {
        assert(c.size() == data.size());

        for (std::size_t i = 0; i < c.size(); i++)
            f(i, c[i], data[i], std::forward<Args_>(args)...);
    }
};

}   // namespace AssemblerLib

#endif  // ASSEMBLERLIB_SERIALEXECUTOR_H_H
