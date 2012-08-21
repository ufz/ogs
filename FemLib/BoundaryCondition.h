/**
 * \file BoundaryCondition.h
 * 2011/08/30 KR inital implementation
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
	BoundaryCondition(const CBoundaryCondition &bc, const std::string &geometry_name);
	BoundaryCondition(const FEMCondition &cond)
		: FEMCondition(cond, FEMCondition::BOUNDARY_CONDITION) {};
	~BoundaryCondition() {}

	size_t getTimType() const {return _tim_type; }
	void setTimType(size_t value) { _tim_type = value; }

private:
	size_t _tim_type;
};

#endif //BOUNDARYCONDITION_H
