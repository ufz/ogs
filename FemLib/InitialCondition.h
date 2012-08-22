/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.net/LICENSE.txt
 *
 * \file InitialCondition.h
 *
 * Created on 2011-08-30 by Karsten Rink
 *
 */

#ifndef INITIALCONDITION_H
#define INITIALCONDITION_H

#include "FEMCondition.h"

/**
 * \brief Adapter class for handling boundary conditions in the user Interface
 * \sa FEMCondition
 */
class InitialCondition : public FEMCondition
{
public:
	InitialCondition(const std::string &geometry_name)
		: FEMCondition(geometry_name, FEMCondition::INITIAL_CONDITION) {};
	//InitialCondition(const CInitialCondition &ic, const std::string &geometry_name);
	InitialCondition(const FEMCondition &cond)
		: FEMCondition(cond, FEMCondition::INITIAL_CONDITION) {};
	~InitialCondition() {}
};

#endif //INITIALCONDITION_H
