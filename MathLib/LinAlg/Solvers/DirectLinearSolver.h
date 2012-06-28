/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.net/LICENSE.txt
 *
 *
 * \file DirectLinearSolver.h
 *
 * Created on 2011-01-07 by Thomas Fischer
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
