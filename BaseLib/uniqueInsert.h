/**
 * \file
 * \author Thomas Fischer
 * \date   2011-02-23
 * \brief  Definition of the uniquePushBack function.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef UNIQUELISTINSERT_H_
#define UNIQUELISTINSERT_H_

#include <algorithm>

namespace BaseLib {

template<typename Container>
void uniquePushBack(Container& container, typename Container::value_type const& element)
{
    if (std::find(container.begin(), container.end(), element) == container.end())
        container.push_back(element);
}

} // end namespace BaseLib

#endif /* UNIQUELISTINSERT_H_ */
