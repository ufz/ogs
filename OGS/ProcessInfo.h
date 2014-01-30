/**
 * \file
 * \author Thomas Fischer
 * \date   2010-09-02
 * \brief  Definition of the ProcessInfo class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESSINFO_H_
#define PROCESSINFO_H_

// FEM
#include "FEMEnums.h"

/**
 * \brief Class ProcessInfo stores the process type,
 * an value for the primary variable of the process and
 * a pointer to the process object.
 */
class ProcessInfo
{
public:
	/**
	 * Default constructor, initializes pcs_type with ProcessType::INVALID_PROCESS,
	 * pcs_pv with PrimaryVariable::INVALID_PV and
	 * the pointer to the process with NULL. The user should set the values
	 * with the appropriate set methods.
	 */
	ProcessInfo();

	/**
	 * constructor initializing all attributes of the object with the given values
	 * @param pcs_type process type (\sa enum ProcessType)
	 * @param pcs_pv type of primary variable (\sa enum PrimaryVariable)
	 * @param pcs a pointer to the process
	 */
	ProcessInfo (FiniteElement::ProcessType pcs_type, FiniteElement::PrimaryVariable pcs_pv/* TODO6 , CRFProcess* pcs*/);

	/**
	 * Sets the process type.
	 * @param pcs_type the process type, for valid values see enum ProcessType
	 */
	void setProcessType (FiniteElement::ProcessType pcs_type);

	/**
	 * Sets the value for the primary variable
	 * @param pcs_pv value for primary variable, possible values are documented in enum PrimaryVariable
	 */
	void setProcessPrimaryVariable (FiniteElement::PrimaryVariable pcs_pv);

	/**
	 * Get the process type.
	 * @return the process type
	 */
	FiniteElement::ProcessType getProcessType () const;

	/**
	 * Get the primary variable of the process.
	 * @return the primary variable of the process
	 */
	FiniteElement::PrimaryVariable getProcessPrimaryVariable () const;

	virtual ~ProcessInfo();

protected:
	/**
	 * process type, see enum ProcessType for valid values
	 */
	FiniteElement::ProcessType _pcs_type;
	/**
	 * the primary variable of the process, see enum PrimaryVariable for valid values
	 */
	FiniteElement::PrimaryVariable _pcs_pv;
};
#endif /* PROCESSINFO_H_ */
