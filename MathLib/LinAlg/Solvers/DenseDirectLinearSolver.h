/**
 * \file DenseDirectLinearSolver.h
 *
 * Created on 2011-01-07 by Thomas Fischer
 */

#ifndef DENSEDIRECTLINEARSOLVER_H_
#define DENSEDIRECTLINEARSOLVER_H_

#include "DirectLinearSolver.h"

namespace MathLib {

class DenseDirectLinearSolver: public MathLib::DirectLinearSolver {
public:
	DenseDirectLinearSolver() {};
	virtual ~DenseDirectLinearSolver() {};
};

}

#endif /* DENSEDIRECTLINEARSOLVER_H_ */
