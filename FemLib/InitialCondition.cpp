/**
 * \file
 * \author Karsten Rink
 * \date   2011-08-30
 * \brief  Implementation of the InitialCondition class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
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
