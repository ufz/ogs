/**
 * \file   IDiscreteVectorBase.h
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

#ifndef IDISCRETEVECTORBASE_H_
#define IDISCRETEVECTORBASE_H_

#include <vector>

#include "DiscreteObjectWithID.h"

namespace DiscreteLib
{

/**
 * \brief Base class for any data types of discrete data
 */
class IDiscreteVectorBase : public DiscreteObjectWithID
{
public:
    virtual ~IDiscreteVectorBase() {};
};

} // end

#endif //IDISCRETEVECTORBASE_H_
