/**
 * \file
 * \author Karsten Rink
 * \date   2011-08-30
 * \brief  Definition of the BoundaryCondition class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef BOUNDARYCONDITION_H
#define BOUNDARYCONDITION_H

#include "FEMCondition.h"

/**
 * \brief Adapter class for handling boundary conditions in the user Interface
 * \sa FEMCondition
 */
class BoundaryCondition : public FEMCondition
{
public:
	BoundaryCondition(const std::string &geometry_name)
		: FEMCondition(geometry_name, FEMCondition::BOUNDARY_CONDITION), _tim_type(0) {};
	//BoundaryCondition(const CBoundaryCondition &bc, const std::string &geometry_name);
	BoundaryCondition(const FEMCondition &cond)
		: FEMCondition(cond, FEMCondition::BOUNDARY_CONDITION) {};
	~BoundaryCondition() {}

	std::size_t getTimType() const {return _tim_type; }
	void setTimType(std::size_t value) { _tim_type = value; }

private:
	std::size_t _tim_type;
};

#endif //BOUNDARYCONDITION_H
