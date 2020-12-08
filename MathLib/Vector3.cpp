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

double scalarTriple(MathLib::Vector3 const& u, MathLib::Vector3 const& v,
                    MathLib::Vector3 const& w)
{
    auto const pu =
        Eigen::Map<Eigen::Vector3d>(const_cast<double*>(u.getCoords()));
    auto const pv =
        Eigen::Map<Eigen::Vector3d>(const_cast<double*>(v.getCoords()));
    auto const pw =
        Eigen::Map<Eigen::Vector3d>(const_cast<double*>(w.getCoords()));
    return pu.cross(pv).dot(pw);
}

//double scalarTriple(EigenLib::Vector3d const& u, EigenLib::Vector3d const& v,
//                    EigenLib::Vector3d const& w)
//{
//    return (u.cross(v).dot(w))(0,0);
//}

}  // end namespace MathLib
