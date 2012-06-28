/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file IterativeLinearSolver.h
 *
 * Created on 2011-01-07 by Thomas Fischer
 */

#ifndef ITERATIVELINEARSOLVER_H_
#define ITERATIVELINEARSOLVER_H_

#include <LinearSolver.h>

namespace MathLib {

class IterativeLinearSolver: public MathLib::LinearSolver {
public:
	IterativeLinearSolver() {};
	virtual ~IterativeLinearSolver() {};
};

}

#endif /* ITERATIVELINEARSOLVER_H_ */
