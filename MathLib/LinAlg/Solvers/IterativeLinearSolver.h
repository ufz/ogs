/*
 * IterativeLinearSolver.h
 *
 *  Created on: Jan 7, 2011
 *      Author: TF
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
