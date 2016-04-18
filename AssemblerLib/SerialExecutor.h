/**
 * \author Norihiro Watanabe
 * \date   2013-04-16
 *
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
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
    /// Executes \c f for each element from the input container.
    ///
    /// The elements of \c c are dereferenced before being passed to \c f.
    /// Return values of the function call are ignored.
    ///
    /// \note
    /// This function is intended to be used if \c f is a plain function,
    /// a static class method, a lambda or an object having a call operator
    /// defined. For non-static member functions you are encouraged to use
    /// executeMemberDereferenced(), if applicable.
    ///
    /// \tparam F    type of the callback function \c f.
    /// \tparam C    input container type.
    /// \tparam Args types additional arguments passed to \c f.
    ///
    /// \param f    a function that accepts a a reference to container's elements,
    ///             an index as arguments (and possibly further arguments).
    /// \param c    a container supporting access over operator[].
    ///             The elements of \c c must have pointer semantics, i.e.,
    ///             support dereferencing via unary operator*().
    /// \param args additional arguments passed to \c f.
    ///
    template <typename F, typename C, typename ...Args_>
    static
    void
#if defined(_MSC_VER) && (_MSC_VER >= 1700)
    executeDereferenced(F& f, C const& c, Args_&&... args)
#else
    executeDereferenced(F const& f, C const& c, Args_&&... args)
#endif
    {
        for (std::size_t i = 0; i < c.size(); i++)
            f(i, *c[i], std::forward<Args_>(args)...);
    }

    /// Executes the given \c method of the given \c object for each element
    /// from the input \c container.
    ///
    /// This method is very similar to executeDereferenced().
    ///
    /// \see executeDereferenced()
    template <typename Container, typename Object, typename... MethodArgs, typename... Args>
    static void
    executeMemberDereferenced(
            Object& object,
            void (Object::*method)(MethodArgs... method_args) const,
            Container const& container, Args&&... args)
    {
        for (std::size_t i = 0; i < container.size(); i++) {
            (object.*method)(i, *container[i], std::forward<Args>(args)...);
        }
    }

    /// \override
    template <typename Container, typename Object, typename... MethodArgs, typename... Args>
    static void
    executeMemberDereferenced(
            Object& object,
            void (Object::*method)(MethodArgs... method_args),
            Container const& container, Args&&... args)
    {
        for (std::size_t i = 0; i < container.size(); i++) {
            (object.*method)(i, *container[i], std::forward<Args>(args)...);
        }
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
    /// \param args additional arguments passed to \c f
    template <typename F, typename C, typename Data, typename ...Args_>
    static
    void
#if defined(_MSC_VER) && (_MSC_VER >= 1700)
    transform(F& f, C const& c, Data& data, Args_&&... args)
#else
    transform(F const& f, C const& c, Data& data, Args_&&... args)
#endif
    {
        assert(c.size() == data.size());

        for (std::size_t i = 0; i < c.size(); i++)
            f(i, c[i], data[i], std::forward<Args_>(args)...);
    }
};

}   // namespace AssemblerLib

#endif  // ASSEMBLERLIB_SERIALEXECUTOR_H_H
