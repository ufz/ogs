/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "Functional.h"

namespace BaseLib
{
namespace detail
{
#define DEFINE_INDEXEDPLACEHOLDER_MEMBER(INDEX, INDEX_P_1) \
    const decltype(std::placeholders::_##INDEX_P_1)        \
        IndexedPlacedPlaceholder<(INDEX)>::value =         \
            std::placeholders::_##INDEX_P_1

DEFINE_INDEXEDPLACEHOLDER_MEMBER(0, 1);
DEFINE_INDEXEDPLACEHOLDER_MEMBER(1, 2);
DEFINE_INDEXEDPLACEHOLDER_MEMBER(2, 3);
DEFINE_INDEXEDPLACEHOLDER_MEMBER(3, 4);
DEFINE_INDEXEDPLACEHOLDER_MEMBER(4, 5);
DEFINE_INDEXEDPLACEHOLDER_MEMBER(5, 6);
DEFINE_INDEXEDPLACEHOLDER_MEMBER(6, 7);
DEFINE_INDEXEDPLACEHOLDER_MEMBER(7, 8);
DEFINE_INDEXEDPLACEHOLDER_MEMBER(8, 9);
DEFINE_INDEXEDPLACEHOLDER_MEMBER(9, 10);
}

}  // namespace BaseLib
