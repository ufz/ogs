/**
 * \file DistributionInfo.cpp
 *
 * Created on 2010-09-28 by Thomas Fischer
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
