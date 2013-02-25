/**
 * \file
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

#include "IDiscreteSystem.h"

#include "IDiscreteLinearSystem.h"
#include "IDiscreteVectorBase.h"

namespace DiscreteLib
{

void IDiscreteSystem::addLinearSystem(IDiscreteLinearSystem *eqs)
{
    _resource.addLinearSystem(eqs);
}

void IDiscreteSystem::deleteLinearSystem(IDiscreteLinearSystem* eqs)
{
    if (eqs!=nullptr) {
        _resource.eraseLinearSystem(eqs);
        delete eqs;
    }
}

void IDiscreteSystem::addVector(IDiscreteVectorBase *vec)
{
    _resource.addVector(vec);
}

void IDiscreteSystem::deleteVector(IDiscreteVectorBase* v)
{
    if (v!=nullptr) {
        _resource.eraseVector(v);
        delete v;
    }
};

} //DiscreteLib

