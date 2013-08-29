/**
 * \file
 * \author Karsten Rink
 * \date   2011-08-30
 * \brief  Implementation of the BoundaryCondition class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */
/*
#include "BoundaryCondition.h"

#include "rf_bc_new.h"

BoundaryCondition::BoundaryCondition(const CBoundaryCondition &bc, const std::string &geometry_name)
	: FEMCondition(geometry_name, bc.getProcessType(), bc.getProcessPrimaryVariable(),
	               bc.getGeomType(), bc.getGeoName(),
	               bc.getProcessDistributionType(), FEMCondition::BOUNDARY_CONDITION)
{
	if (this->getProcessDistributionType() == FiniteElement::CONSTANT ||
	    this->getProcessDistributionType() == FiniteElement::CONSTANT_NEUMANN)
		this->setConstantDisValue(bc.getGeoNodeValue());
	else if (this->getProcessDistributionType() == FiniteElement::LINEAR ||
	         this->getProcessDistributionType() == FiniteElement::LINEAR_NEUMANN)
	{
		const std::vector<int> bc_nodes(bc.getPointsWithDistribedBC());
		std::vector<size_t> dis_nodes(bc_nodes.size());
		for (size_t i=0; i<dis_nodes.size(); i++) dis_nodes[i] = static_cast<size_t>(bc_nodes[i]);
		this->setDisValues(dis_nodes, bc.getDistribedBC());
	}
}
*/
