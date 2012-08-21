/**
 * \file InitialCondition.cpp
 * 2011/08/30 KR inital implementation
 *
 */

#include "InitialCondition.h"
#include "rf_ic_new.h"

InitialCondition::InitialCondition(const CInitialCondition &ic, const std::string &geometry_name)
	: FEMCondition(geometry_name, ic.getProcessType(), ic.getProcessPrimaryVariable(),
	               ic.getGeoType(),
	               (ic.getGeoType() == GEOLIB::GEODOMAIN) ? "Domain" : ic.getGeoName(),
	               ic.getProcessDistributionType(), FEMCondition::INITIAL_CONDITION)
{
	if (this->getProcessDistributionType() == FiniteElement::CONSTANT)
		this->setConstantDisValue(ic.getGeoNodeValue());
}
