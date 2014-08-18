/**
 * \date   2014-08-04
 * \brief  Implementation of OpenGeoSys simulation application
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <cstdlib>

// ThirdParty/logog
#include "logog/include/logog.hpp"

// ThirdParty/tclap
#include "ThirdParty/tclap/CmdLine.h"

// BaseLib
#include "BaseLib/BuildInfo.h"
#include "BaseLib/LogogSimpleFormatter.h"

// OGS

#include "Applications/ApplicationsLib/ProjectData.h"

int main(int argc, char *argv[])
{
	// logog
	LOGOG_INITIALIZE();
	BaseLib::LogogSimpleFormatter *fmt(new BaseLib::LogogSimpleFormatter);
	logog::Cout *logog_cout(new logog::Cout);
	logog_cout->SetFormatter(*fmt);

	// Parse CLI arguments.
	TCLAP::CmdLine cmd("OpenGeoSys-6 software (" +
			BaseLib::BuildInfo::git_version_sha1_short + ")\n" +
			"Copyright (c) 2014, OpenGeoSys Community " +
			"(http://www.opengeosys.org) " +
			"Distributed under a Modified BSD License. " +
			"See accompanying file LICENSE.txt or " +
			"http://www.opengeosys.org/project/license",
		' ',
		BaseLib::BuildInfo::git_version_sha1);

	TCLAP::ValueArg<std::string> project_arg(
		"p", // option tag
		"project", // long option tag
		"file name of the project", // description
		true, // required
		"",
		"string"); // type of information
	cmd.add(project_arg);
	cmd.parse(argc, argv);

	ProjectData project_data;

	delete fmt;
	delete logog_cout;
	LOGOG_SHUTDOWN();

	return 0;
}
