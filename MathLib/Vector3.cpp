/**
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "Vector3.h"

namespace MathLib
{

double scalarTriple(MathLib::Vector3 const& u, MathLib::Vector3 const& v,
                    MathLib::Vector3 const& w)
{
    MathLib::Vector3 const cross(MathLib::crossProduct(u, v));
    return MathLib::scalarProduct(cross,w);
}

}  // end namespace MathLib
