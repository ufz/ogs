/**
 * \file
 * \author Thomas Fischer
 * \date   2010-09-02
 * \brief  Implementation of the ProcessInfo class.
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <ProcessInfo.h>

ProcessInfo::ProcessInfo() :
	_pcs_type (FiniteElement::INVALID_PROCESS), _pcs_pv (FiniteElement::INVALID_PV)
{}

ProcessInfo::ProcessInfo (FiniteElement::ProcessType pcs_type, FiniteElement::PrimaryVariable pcs_pv) :
	_pcs_type (pcs_type), _pcs_pv (pcs_pv)
{}

void ProcessInfo::setProcessType (FiniteElement::ProcessType pcs_type)
{
	_pcs_type = pcs_type;
}

void ProcessInfo::setProcessPrimaryVariable (FiniteElement::PrimaryVariable pcs_pv)
{
	_pcs_pv = pcs_pv;
}

FiniteElement::ProcessType ProcessInfo::getProcessType () const
{
	return _pcs_type;
}

FiniteElement::PrimaryVariable ProcessInfo::getProcessPrimaryVariable () const
{
	return _pcs_pv;
}

ProcessInfo::~ProcessInfo()
{}
