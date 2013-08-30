/**
 * \author Norihiro Watanabe
 * \date   2013-08-13
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */


#include "GaussLegendre.h"

namespace MathLib
{

std::pair<double, double> GaussLegendre::getPoint(unsigned n_sample_points, unsigned point_id)
{
    switch (n_sample_points)
    {
    case 1:
        return std::make_pair(0.0, 2.0);
    case 2:
        switch (point_id)
        {
        case 0:
            return std::make_pair(0.577350269189626, 1.0);
        case 1:
            return std::make_pair(-0.577350269189626, 1.0);
        }
        break;
    case 3:
        switch (point_id)
        {
        case 0:
            return std::make_pair(0.774596669241483, 0.555555555555556);
        case 1:
            return std::make_pair(0.0, 0.888888888888889);
        case 2:
            return std::make_pair(-0.774596669241483, 0.555555555555556);
        }
        break;
    case 4:
        switch (point_id)
        {
        case 0:
            return std::make_pair(0.861136311594053, 0.347854845137454);
        case 1:
            return std::make_pair(0.339981043584856, 0.652145154862546);
        case 2:
            return std::make_pair(-0.339981043584856, 0.652145154862546);
        case 3:
            return std::make_pair(-0.861136311594053, 0.347854845137454);
        }
        break;
    default:
        return std::make_pair(0.0, 0.0);
    }
}


} //namespace
