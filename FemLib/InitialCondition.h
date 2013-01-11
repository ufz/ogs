/**
 * \file
 * \author Karsten Rink
 * \date   2011-08-30
 * \brief  Definition of the InitialCondition class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
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
