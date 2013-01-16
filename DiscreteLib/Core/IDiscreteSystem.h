/**
 * \file   IDiscreteSystem.h
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

#ifndef IDISCRETESYSTEM_H_
#define IDISCRETESYSTEM_H_

#include "DiscreteResourceManager.h"

namespace MeshLib
{
class Mesh;
}

namespace DiscreteLib
{

// forward declaration
class IDiscreteLinearSystem;
class IDiscreteVectorBase;

/**
 * \brief Interface of discrete systems
 *
 *  Discrete system class contains the followings
 *  - discrete space (i.e. mesh)
 *  - discrete data (i.e. vector)
 *  - linear systems which are used to calculate discrete data
 */
class IDiscreteSystem
{
public:
    virtual ~IDiscreteSystem() {};
    
    /// return a reference to a mesh object owned by this system
    virtual const MeshLib::Mesh& getMesh() const = 0;

protected:
    /// add an equation object to this system
    void addLinearSystem(IDiscreteLinearSystem *eqs);

    /// delete an equation object from this system
    void deleteLinearSystem(IDiscreteLinearSystem* eqs);

    /// add a vector object to this system
    void addVector(IDiscreteVectorBase *vec);

    /// delete a vector object from this system
    void deleteVector(IDiscreteVectorBase* vec);

    DiscreteResourceManager _resource;
};

} //end

#endif //IDISCRETESYSTEM_H_
