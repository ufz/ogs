/**
 * @file BoostXmlCndInterface.cpp
 * @author git blame BoostXmlCndInterface.cpp
 * @date Oct 14, 2013
 * @brief 
 *
 * @copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */
#include "BoostXmlCndInterface.h"

#include <fstream>

#include "logog/include/logog.hpp"

#include "BoundaryCondition.h"

namespace FileIO
{

BoostXmlCndInterface::BoostXmlCndInterface(ProjectData* project_data) :
		_project_data(project_data)
{}

bool BoostXmlCndInterface::readFile(const std::string &fname)
{
	return false;
}

} // end namespace FileIO
