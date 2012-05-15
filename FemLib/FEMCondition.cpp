/**
 * \file FEMCondition.cpp
 * 25/11/2010 KR inital implementation
 *
 */

#include "FEMCondition.h"
#include "ProcessInfo.h"

#include "GEOObjects.h" //for SourceTerm
#include "GridAdapter.h"

FEMCondition::FEMCondition(const std::string &geometry_name, CondType t)
	: _type(t), _geoName("[unspecified]"), _associated_geometry(geometry_name)
{
	this->setProcessType(FiniteElement::INVALID_PROCESS);
	this->setProcessPrimaryVariable(FiniteElement::INVALID_PV);
	this->setGeoType(GeoLib::INVALID);
	this->setProcessDistributionType(FiniteElement::INVALID_DIS_TYPE);
}

FEMCondition::FEMCondition(const std::string &geometry_name, FiniteElement::ProcessType pt,
				FiniteElement::PrimaryVariable pv, GeoLib::GEOTYPE gt, const std::string &gn,
				FiniteElement::DistributionType dt, CondType ct)
	: ProcessInfo(pt, pv, NULL),  GeoInfo(gt, NULL), DistributionInfo(dt), _type(ct),
	  _geoName(gn), _associated_geometry(geometry_name)
{
}

FEMCondition::FEMCondition(const FEMCondition &cond, CondType t)
	: ProcessInfo(cond.getProcessType(), cond.getProcessPrimaryVariable(), NULL),
	  GeoInfo(cond.getGeoType(), cond.getGeoObj()),
	  DistributionInfo(cond.getProcessDistributionType()),
	  _type(t),
	  _geoName(cond.getGeoName()),
	  _disNodes(cond.getDisNodes()),
	  _disValues(cond.getDisValues()),
	  _associated_geometry(cond.getAssociatedGeometryName())
{
}

std::string FEMCondition::condTypeToString(CondType type)
{
	if (type == FEMCondition::BOUNDARY_CONDITION)
		return "Boundary Conditions";
	else if (type == FEMCondition::INITIAL_CONDITION)
		return "Initial Conditions";
	else if (type == FEMCondition::SOURCE_TERM)
		return "Source Terms";
	else
		return "Unspecified";
}

void FEMCondition::setDisValues(const std::vector< std::pair<size_t, double> > &dis_values)
{
	std::vector<size_t> nodes;
	std::vector<double> values;
	for (size_t i = 0; i < dis_values.size(); i++)
	{
		nodes.push_back(dis_values[i].first);
		values.push_back(dis_values[i].second);
	}
	this->_disNodes=nodes;
	this->_disValues=values;
}
