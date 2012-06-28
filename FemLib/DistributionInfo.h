/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www./**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.net/LICENSE.txt
 *
 *
opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www./**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.net/LICENSE.txt
 *
 *
opengeosys.com/LICENSE.txt
 *
 *
 * \file DistributionInfo.h
 *
 * Created on 2010-09-28 by Thomas Fischer
 */

#ifndef DISTRIBUTIONINFO_H_
#define DISTRIBUTIONINFO_H_

// FEM
#include "FEMEnums.h"

class DistributionInfo
{
public:
	DistributionInfo(FiniteElement::DistributionType dt = FiniteElement::INVALID_DIS_TYPE);
	virtual ~DistributionInfo();

	/**
	 * Sets the value for the distribution type
	 * @param dis_type value for primary variable, possible values are documented in enum PrimaryVariable
	 */
	void setProcessDistributionType (FiniteElement::DistributionType dis_type);

	/**
	 * Get the distribution type of the process.
	 * @return the distribution type of the process
	 */
	FiniteElement::DistributionType getProcessDistributionType () const;

private:
	/**
	 * the distribution type of the process, see enum DistributionType for valid values
	 */
	FiniteElement::DistributionType _dis_type;
};

#endif                                            /* DISTRIBUTIONINFO_H_ */
