/**
 * \file
 * \author Thomas Fischer
 * \date   2011-01-07
 * \brief  Definition of the DirectLinearSolver class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef DIRECTLINEARSOLVER_H_
#define DIRECTLINEARSOLVER_H_

#include "LinearSolver.h"

namespace MathLib {

class DirectLinearSolver : public MathLib::LinearSolver
{
public:
	DirectLinearSolver() {};
	virtual ~DirectLinearSolver() {};
};

}

#endif /* DIRECTLINEARSOLVER_H_ */
