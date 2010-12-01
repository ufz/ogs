/**
 * \file FEMCondition.h
 * 25/11/2010 KR inital implementation
 *
 */

#ifndef FEMCONDITION_H
#define FEMCONDITION_H

#include "GeoInfo.h"
#include "ProcessInfo.h"
#include "DistributionInfo.h"



/** 
 * \brief Adapter class for handling Initial Conditions in the user Interface
 */
class FEMCondition : public GeoInfo, public ProcessInfo, public DistributionInfo
{
public:
	FEMCondition() {};
	~FEMCondition() {};
};

/** 
 * \brief Adapter class for handling Boundary Conditions in the user Interface
 */
class BoundaryCondition : public FEMCondition
{
public:
	BoundaryCondition() : _tim_type(0) {};
	~BoundaryCondition() {};

	size_t getTimType() {return _tim_type; };
	void setTimType(size_t value) { _tim_type = value; };


private:
	size_t _tim_type;
};

/** 
 * \brief Adapter class for handling Source Terms in the user Interface
 */
class InitialCondition : public FEMCondition
{
public:
	InitialCondition() {};
	~InitialCondition() {};
};

/** 
 * \brief Adapter class for handling FEM Conditions in the user Interface
 */
class SourceTerm : public FEMCondition
{
public:
	SourceTerm() : _tim_type(0) {};
	~SourceTerm() {};

	size_t getTimType() {return _tim_type; };
	void setTimType(size_t value) { _tim_type = value; };

private:
	size_t _tim_type;
};


#endif //FEMCONDITION_H
