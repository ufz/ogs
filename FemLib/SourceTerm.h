/**
 * \file
 * \author Karsten Rink
 * \date   2011-08-30
 * \brief  Definition of the SourceTerm class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef SOURCETERM_H
#define SOURCETERM_H

#include "FEMCondition.h"

/**
 * \brief Adapter class for handling boundary conditions in the user Interface
 * \sa FEMCondition
 */
class SourceTerm : public FEMCondition
{
public:
	SourceTerm(const std::string &geometry_name)
		: FEMCondition(geometry_name, FEMCondition::SOURCE_TERM), _tim_type(0) {}
	//SourceTerm(const CSourceTerm &st, const std::string &geometry_name);
	SourceTerm(const FEMCondition &cond)
		: FEMCondition(cond, FEMCondition::SOURCE_TERM) {};
	~SourceTerm() {}

	std::size_t getTimType() const {return _tim_type; }
	void setTimType(std::size_t value) { _tim_type = value; }

	// Legacy function (only required for ascii st-files): reads values for 'direct' source terms
	static void getDirectNodeValues(const std::string &filename,
	                                std::vector< std::pair<std::size_t, double> > &nodes_values);

private:
	std::size_t _tim_type;
};

#endif //SOURCETERM_H
