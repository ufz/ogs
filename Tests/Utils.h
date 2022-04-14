/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

// Used to distinguish function overloads
template <typename T>
struct Type
{
};

// helper type for converting lists of types
template <typename InputList, template <typename...> typename NewListType>
struct ConvertListType;
template <template <typename...> typename OldListType, typename... Ts,
          template <typename...> typename NewListType>
struct ConvertListType<OldListType<Ts...>, NewListType>
{
    using type = NewListType<Ts...>;
};
template <typename InputList, template <typename...> typename NewListType>
using ConvertListType_t =
    typename ConvertListType<InputList, NewListType>::type;

#if __cpp_concepts >= 201907L && __cpp_lib_concepts >= 202002L
#define OGS_HAVE_CONCEPTS
#define OGS_USE_CONCEPT(name) name
#else
#define OGS_USE_CONCEPT(name) typename
#endif
