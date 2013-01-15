/**
 * \file   IDiscreteSystem.h
 * \author Norihiro Watanabe
 * \date   2012-08-03
 * \brief  Helper macros.
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

#include "MeshLib/Mesh.h"
#include "DiscreteDataContainer.h"

namespace DiscreteLib
{

// forward declaration
class IDiscreteLinearSystem;
class IDiscreteVectorBase;

/**
 * \brief Interface for all kinds of discrete systems
 *
 *  Discrete system class contains the followings
 *  - discrete space (i.e. mesh)
 *  - discrete data (i.e. vector)
 *  - linear equations which are used to calculate discrete data
 */
class IDiscreteSystem
{
public:
    virtual ~IDiscreteSystem() {};
    
    /// return a reference to a mesh object owned by this system
    virtual const MeshLib::Mesh& getMesh() const = 0;

    /// add an equation object to this system
    void addLinearSystem(IDiscreteLinearSystem *eqs);

    /// delete an equation object from this system
    void deleteLinearSystem(IDiscreteLinearSystem* eqs);

    /// add a vector object to this system
    void addVector(IDiscreteVectorBase *vec);

    /// delete a vector object from this system
    void deleteVector(IDiscreteVectorBase* vec);

private:
    DiscreteDataContainer _data;
};

void IDiscreteSystem::addLinearSystem(IDiscreteLinearSystem *eqs)
{
    _data.addLinearSystem(eqs);
}

void IDiscreteSystem::deleteLinearSystem(IDiscreteLinearSystem* eqs)
{
    if (eqs!=nullptr) {
        _data.eraseLinearSystem(eqs);
        delete eqs;
    }
}

void IDiscreteSystem::addVector(IDiscreteVectorBase *vec)
{
    _data.addVector(vec);
}

void IDiscreteSystem::deleteVector(IDiscreteVectorBase* v)
{
    if (v!=nullptr) {
        _data.eraseVector(v);
        delete v;
    }
};
} //end

#endif //IDISCRETESYSTEM_H_
