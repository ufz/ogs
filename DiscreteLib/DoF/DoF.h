/**
 * \file   DoF.h
 * \author Norihiro Watanabe
 * \date   2012-08-03
 * \brief
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef DOF_H_
#define DOF_H_

namespace DiscreteLib
{

/**
 * \brief Degree of Freedom (DoF)
 *
 * A DoF is defined as a unique combination of the followings
 * - Variable ID
 * - Mesh ID
 * - discrete point id (Point ID) in the mesh (e.g. node id, element id)
 */
struct DoF
{
    typedef std::size_t VariableID;
    typedef std::size_t MeshID;
    typedef std::size_t PointID;

    VariableID varID;
    MeshID mshID;
    PointID ptID;

    /**
     * constructor
     * @param var   Variable ID
     * @param msh   Mesh ID
     * @param pt    Discrete point ID, e.g. node id or element id
     */
    DoF(VariableID var, MeshID msh, PointID pt)
    : varID(var), mshID(msh), ptID(pt) {};
};


} //end

#endif //DOF_H_
