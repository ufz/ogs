/**
 * \file
 * \author Thomas Fischer
 * \date   2011-01-25
 * \brief  Definition of the ClosestPair class.
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
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
        pnts_ (pnts), id0_ (id0), id1_ (id1)
    {}

protected:
    std::vector<GeoLib::Point*> const & pnts_;
    std::size_t id0_;
    std::size_t id1_;
};

} // end namespace GeoLib
