/**
 * \file   DiscreteResourceManager.h
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

#ifndef DISCRETERESOURCEMANAGER_H_
#define DISCRETERESOURCEMANAGER_H_


#include <vector>

#include "BaseLib/CodingTools.h"

namespace DiscreteLib
{

//forward declaration
class IDiscreteVectorBase;
class IDiscreteLinearSystem;

/**
 * \brief Data container for any discrete objects
 *
 * This class manages discrete vector and linear system objects.
 */
class DiscreteResourceManager
{
public:
    /**
     *
     */
    DiscreteResourceManager(){};

    /**
     *
     */
    virtual ~DiscreteResourceManager();

    /**
     * add a vector
     * @param v
     * @return the vector index
     */
    std::size_t addVector(IDiscreteVectorBase* v);

    /**
     * erase a vector
     * @param v
     */
    void eraseVector(IDiscreteVectorBase* v);

    /**
     * return the number of vectors
     * @return
     */
    std::size_t getNumberOfVectors() const;

    /**
     * return a vector with the given index
     * @param i
     * @return
     */
    IDiscreteVectorBase* getVector(std::size_t i);

    /**
     * add a linear equation object
     * @param eq
     * @return
     */
    std::size_t addLinearSystem(IDiscreteLinearSystem* eq);

    /**
     * erase a linear equation object
     * @param eq
     */
    void eraseLinearSystem(IDiscreteLinearSystem* eq);

    /**
     * return the number of linear equations
     * @return
     */
    std::size_t getNumberOfLinearSystems() const;

    /**
     * return a linear equation with the given index
     * @param i
     * @return
     */
    IDiscreteLinearSystem* getLinearSystem(std::size_t i);

private:
    DISALLOW_COPY_AND_ASSIGN(DiscreteResourceManager);

private:
    std::vector<IDiscreteLinearSystem*> _vec_linear_sys;
    std::vector<IDiscreteVectorBase*> _vec_vectors;
};

} //end

#endif //DISCRETEDATACONTAINER_H_
