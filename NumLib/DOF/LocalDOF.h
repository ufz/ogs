/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <boost/mp11.hpp>

#include "NumLib/Fem/ShapeMatrixPolicy.h"

namespace NumLib
{
template <typename ShapeFunction, int Dim>
struct Vectorial
{
};

namespace detail
{
template <typename ShapeFunction>
struct NumberOfDofs : std::integral_constant<unsigned, ShapeFunction::NPOINTS>
{
};

template <typename ShapeFunction, int Dim>
struct NumberOfDofs<Vectorial<ShapeFunction, Dim>>
    : std::integral_constant<unsigned, ShapeFunction::NPOINTS * Dim>
{
};

template <typename... Ns_t,
          typename ElementDOFVector,
          int... Offsets,
          int... Sizes>
auto localDOFImpl(ElementDOFVector const& x,
                  boost::mp11::mp_list_c<int, Offsets...>,
                  boost::mp11::mp_list_c<int, Sizes...>)
{
    static_assert(((Sizes > 0) && ...));
    static_assert(((Offsets >= 0) && ...));

    return std::tuple<Eigen::Map<typename MatrixPolicyType::VectorType<
        NumberOfDofs<Ns_t>::value> const>...>{{x.data() + Offsets, Sizes}...};
}
}  // namespace detail

//! Decomposes the passed d.o.f.s of an entire finite element into d.o.f.s
//! belonging to individual unknowns—such as pressure, temperature, or
//! displacement—via the passed shape functions.
//!
//! \attention The order of the passed shape functions indicates the order of
//! the primary unknowns in the element d.o.f.s, and the order of the returned
//! data, i.e., the order is important!
//!
//! \return a std::tuple of mapped d.o.f. data (cf. Eigen::Map)
template <typename... Ns_t, typename ElementDOFVector>
auto localDOF(ElementDOFVector const& x)
{
    using namespace boost::mp11;

    static_assert(sizeof...(Ns_t) > 0);

    using Sizes = mp_list_c<int, detail::NumberOfDofs<Ns_t>::value...>;
    using Offsets =
        mp_push_front<mp_pop_back<mp_partial_sum<Sizes, mp_int<0>, mp_plus>>,
                      mp_int<0>>;

    assert((mp_back<Offsets>::value + mp_back<Sizes>::value == x.size()) &&
           "The passed shape matrices require a different number of "
           "d.o.f.s than present in the passed element d.o.f. vector.");

    return detail::localDOFImpl<Ns_t...>(x, Offsets{}, Sizes{});
}

}  // namespace NumLib
