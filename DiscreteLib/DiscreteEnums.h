/**
 * \file   DiscreteEnums.h
 * \author Norihiro Watanabe
 * \date   2012-08-17
 * \brief  Enum definitions
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef DISCRETEENUMS_H_
#define DISCRETEENUMS_H_

#include <valarray>
#include "MathLib/LinAlg/Dense/Matrix.h"
#include "MathLib/LinAlg/SystemOfLinearEquations/DenseLinearSystem.h"

namespace DiscreteLib
{

/**
 * \brief Discrete System Type
 */
struct DiscreteSystemType
{
    enum type
    {
        Serial
    };
};

typedef MathLib::DenseLinearSystem LocalLinearSystem;
typedef LocalLinearSystem::MatrixType LocalMatrix;
typedef LocalLinearSystem::VectorType LocalVector;


}

#endif //DISCRETEENUMS_H_
