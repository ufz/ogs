/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www./**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.net/LICENSE.txt
 *
 *
opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www./**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.net/LICENSE.txt
 *
 *
opengeosys.com/LICENSE.txt
 *
 *
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
