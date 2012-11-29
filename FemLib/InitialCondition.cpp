/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 *
 * \file InitialCondition.cpp
 *
 * Created on 2011-08-30 by Karsten Rink
 */

/*
#include "InitialCondition.h"
#include "rf_ic_new.h"

InitialCondition::InitialCondition(const CInitialCondition &ic, const std::string &geometry_name)
	: FEMCondition(geometry_name, ic.getProcessType(), ic.getProcessPrimaryVariable(),
	               ic.getGeomType(),
	               (ic.getGeomType() == GeoLib::GEODOMAIN) ? "Domain" : ic.getGeoName(),
	               ic.getProcessDistributionType(), FEMCondition::INITIAL_CONDITION)
{
	if (this->getProcessDistributionType() == FiniteElement::CONSTANT)
		this->setConstantDisValue(ic.getGeoNodeValue());
}
*/
