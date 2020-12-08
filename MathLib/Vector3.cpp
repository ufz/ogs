/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "Vector3.h"

namespace MathLib
{

double scalarTriple(Eigen::Vector3d const& u, Eigen::Vector3d const& v,
                    Eigen::Vector3d const& w)
{
    return u.cross(v).dot(w);
}

}  // end namespace MathLib
