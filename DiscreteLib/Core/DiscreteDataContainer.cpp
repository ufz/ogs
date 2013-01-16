/**
 * \file   DiscreteDataContainer.cpp
 * \author Norihiro Watanabe
 * \date   2012-08-03
 * \brief  Discrete data container
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "DiscreteDataContainer.h"

#include "IDiscreteVector.h"
#include "IDiscreteLinearSystem.h"

namespace DiscreteLib
{

DiscreteDataContainer::~DiscreteDataContainer()
{
    BaseLib::releaseObjectsInStdVector(_vec_vectors);
    BaseLib::releaseObjectsInStdVector(_vec_linear_sys);
}

std::size_t DiscreteDataContainer::addVector(IDiscreteVectorBase* v)
{
    _vec_vectors.push_back(v);
    v->setObjectID(_vec_vectors.size()-1);
    return _vec_vectors.size()-1;
}

void DiscreteDataContainer::eraseVector(IDiscreteVectorBase* v)
{
    const std::size_t i = v->getObjectID();
    if (_vec_vectors.size() > i) {
        _vec_vectors[i] = nullptr;
    }
}

std::size_t DiscreteDataContainer::getNumberOfVectors() const
{
    return _vec_vectors.size();
};

IDiscreteVectorBase* DiscreteDataContainer::getVector(std::size_t i)
{
    return _vec_vectors[i];
}

std::size_t DiscreteDataContainer::addLinearSystem(IDiscreteLinearSystem* eq)
{
    _vec_linear_sys.push_back(eq);
    eq->setObjectID(_vec_linear_sys.size()-1);
    return _vec_linear_sys.size()-1;
}

void DiscreteDataContainer::eraseLinearSystem(IDiscreteLinearSystem* eq)
{
    const std::size_t i = eq->getObjectID();
    if (_vec_linear_sys.size() > i) {
        _vec_linear_sys[i] = nullptr;
    }
}

std::size_t DiscreteDataContainer::getNumberOfLinearSystems() const
{
    return _vec_linear_sys.size();
};

IDiscreteLinearSystem* DiscreteDataContainer::getLinearSystem(std::size_t i)
{
    return _vec_linear_sys[i];
}


} //end
