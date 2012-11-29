/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 *
 * \file SourceTerm.cpp
 *
 * Created on 2011-08-30 by Karsten Rink
 *
 */
/*
#include "SourceTerm.h"
#include "rf_st_new.h"

SourceTerm::SourceTerm(const CSourceTerm &st, const std::string &geometry_name)
	: FEMCondition(geometry_name, st.getProcessType(), st.getProcessPrimaryVariable(),
	               st.getGeomType(), st.getGeoName(),
	               st.getProcessDistributionType(), FEMCondition::SOURCE_TERM)
{
	if (this->getProcessDistributionType() == FiniteElement::CONSTANT ||
	    this->getProcessDistributionType() == FiniteElement::CONSTANT_NEUMANN)
		this->setConstantDisValue(st.getGeoNodeValue());
	else if (this->getProcessDistributionType() == FiniteElement::LINEAR ||
	         this->getProcessDistributionType() == FiniteElement::LINEAR_NEUMANN)
	{
		const std::vector<int> st_nodes(st.getPointsWithDistribedST());
		std::vector<size_t> dis_nodes(st_nodes.size());
		for (size_t i=0; i<dis_nodes.size(); i++) dis_nodes[i] = static_cast<size_t>(st_nodes[i]);
		this->setDisValues(dis_nodes, st.getDistribedST());
	}
	else if (this->getProcessDistributionType() == FiniteElement::DIRECT)
	{
		//this->_direct_file_name = st.fname;
	}
	else
		std::cout << "Error in SourceTerm() - Unknown Process Distribution Type \"" <<
		FiniteElement::convertDisTypeToString(st.getProcessDistributionType()) <<
		"\"..." <<
		std::endl;
}

// Legacy function (only required for ascii st-files): reads values for 'direct' source terms
void SourceTerm::getDirectNodeValues(const std::string &filename,
                                     std::vector< std::pair<size_t, double> > &node_values)
{
	std::ifstream in(filename.c_str());
	if (!in.is_open())
	{
		std::cout << "Error in getNodeValues() - Could not find file for direct node values..." << std::endl;
		return;
	}

	std::stringstream str_in;
	std::string line("");
	size_t idx(0);
	double val(0);

	while ( getline(in, line) )
	{
		if (line.find("#STOP") != std::string::npos)
			return;
		str_in << line;
		str_in >> idx >> val;
		node_values.push_back(std::pair<size_t, double>(idx, val));
		str_in.clear();
	}
}
*/