/**
 * \file
 * \author Thomas Fischer
 * \date   2010-09-28
 * \brief  Implementation of the DistributionInfo class.
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "DistributionInfo.h"

DistributionInfo::DistributionInfo(FiniteElement::DistributionType dt) :
	_dis_type (dt)
{}

DistributionInfo::~DistributionInfo()
{}

void DistributionInfo::setProcessDistributionType (FiniteElement::DistributionType dis_type)

{
	_dis_type = dis_type;
}

FiniteElement::DistributionType DistributionInfo::getProcessDistributionType () const
{
	return _dis_type;
}
