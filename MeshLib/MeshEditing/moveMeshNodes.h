/**
 * @file moveMeshNodes.h
 * @author Thomas Fischer
 * @date Feb 03, 2014
 * @brief Functionality to move mesh nodes using a given displacement vec.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */
#pragma once

#include "MeshLib/Node.h"

namespace MeshLib
{

/**
    * Function that moves mesh nodes.
    *
    * The function iterates over all mesh nodes between the
    * begin and the end iterator and moves them using the
    * given displacement.
    * @param begin begin iterator
    * @param end end iterator
    * @param displacement the displacement to use
*/
template <typename Iterator>
void moveMeshNodes(
    Iterator begin,
    Iterator end,
    MeshLib::Node const& displacement)
{
    std::for_each(begin, end, [&displacement](MeshLib::Node* node)
        {
            (*node)[0] += displacement[0];
            (*node)[1] += displacement[1];
            (*node)[2] += displacement[2];
        }
    );
};

} // end namespace MeshLib
