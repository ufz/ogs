/**
 * \file
 * \author Thomas Fischer
 * \date   2011-01-07
 * \brief  Definition of the DenseDirectLinearSolver class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
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
