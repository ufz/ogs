/**
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
