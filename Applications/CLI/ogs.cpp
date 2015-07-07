/**
 * \date   2014-08-04
 * \brief  Implementation of OpenGeoSys simulation application
 *
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

// ThirdParty/logog
#include "logog/include/logog.hpp"

// ThirdParty/tclap
#include "ThirdParty/tclap/CmdLine.h"

// BaseLib
#include "BaseLib/BuildInfo.h"
#include "BaseLib/FileTools.h"
#include "BaseLib/LogogSimpleFormatter.h"

#include "Applications/ApplicationsLib/ProjectData.h"

#include "NumericsConfig.h"


int main(int argc, char *argv[])
{

	using ConfigTree = boost::property_tree::ptree;

	// logog
	LOGOG_INITIALIZE();
	BaseLib::LogogSimpleFormatter *fmt(new BaseLib::LogogSimpleFormatter);
	logog::Cout *logog_cout(new logog::Cout);
	logog_cout->SetFormatter(*fmt);

	// Parse CLI arguments.
	TCLAP::CmdLine cmd("OpenGeoSys-6 software.\n"
			"Copyright (c) 2012-2015, OpenGeoSys Community "
			"(http://www.opengeosys.org) "
			"Distributed under a Modified BSD License. "
			"See accompanying file LICENSE.txt or "
			"http://www.opengeosys.org/project/license",
		' ',
		BaseLib::BuildInfo::git_describe);

	TCLAP::UnlabeledValueArg<std::string> project_arg(
		"project-file",
		"Path to the ogs6 project file.",
		true,
		"",
		"PROJECT FILE");

	cmd.add(project_arg);
	cmd.parse(argc, argv);

#ifdef USE_LIS
	lis_initialize(&argc, &argv);
#endif

	// Project's configuration
	ConfigTree project_config;

	read_xml(project_arg.getValue(), project_config,
			boost::property_tree::xml_parser::no_comments);
	DBUG("Project configuration from file \'%s\' read.",
		project_arg.getValue().c_str());


	project_config = project_config.get_child("OpenGeoSysProject");

	ProjectData project(project_config,
			BaseLib::extractPath(project_arg.getValue()));

	// Create processes.
	project.buildProcesses<GlobalSetupType>();

	INFO("Initialize processes.");
	for (auto p_it = project.processesBegin(); p_it != project.processesEnd(); ++p_it)
	{
		(*p_it)->initialize();
	}

	std::string const output_file_name(project.getOutputFilePrefix() + ".vtu");

	INFO("Solve processes.");
	for (auto p_it = project.processesBegin(); p_it != project.processesEnd(); ++p_it)
	{
		(*p_it)->solve();
		(*p_it)->post(output_file_name);
	}

	delete fmt;
	delete logog_cout;
	LOGOG_SHUTDOWN();

#ifdef USE_LIS
	lis_finalize();
#endif

	return 0;
}
