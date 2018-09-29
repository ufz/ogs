/**
 * \file
 * \author Thomas Fischer
 * \date   2011-01-25
 * \brief  Definition of the ClosestPair class.
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

// STL
#include <vector>

// GeoLib
#include "Point.h"

namespace GeoLib {

class ClosestPair
{
public:
    ClosestPair (std::vector<GeoLib::Point*> const & pnts, std::size_t id0, std::size_t id1) :
        _pnts (pnts), _id0 (id0), _id1 (id1)
    {}

protected:
    std::vector<GeoLib::Point*> const & _pnts;
    std::size_t _id0;
    std::size_t _id1;
};

} // end namespace GeoLib
