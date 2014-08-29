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

#include "FileIO/readMeshFromFile.h"
#include "FileIO/XmlIO/Boost/BoostXmlGmlInterface.h"

#include "GeoLib/GEOObjects.h"

#include "MeshLib/Mesh.h"

#include "ProcessLib/ProcessVariable.h"
#include "ProcessLib/Process.h"

#include "Applications/ApplicationsLib/ProjectData.h"



int main(int argc, char *argv[])
{

	using ConfigTree = boost::property_tree::ptree;

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

	TCLAP::UnlabeledValueArg<std::string> project_arg(
		"project-file",
		"Path to the ogs6 project file.",
		true,
		"",
		"PROJECT FILE");

	cmd.add(project_arg);
	cmd.parse(argc, argv);

	// Project's configuration
	ConfigTree project_config;

	read_xml(project_arg.getValue(), project_config, false);
	DBUG("Project configuration from file \'%s\' read.",
		project_arg.getValue().c_str());


	project_config = project_config.get_child("OpenGeoSysProject");

	// Mesh
	MeshLib::Mesh const* mesh;
	{
		std::string const mesh_file =
			BaseLib::copyPathToFileName(
					project_config.get<std::string>("mesh"),
					project_arg.getValue());
		mesh = FileIO::readMeshFromFile(mesh_file);
	}

	// Geometry
	GeoLib::GEOObjects geometries;
	{
		std::string const geometry_file =
			BaseLib::copyPathToFileName(
					project_config.get<std::string>("geometry"),
					project_arg.getValue());
		DBUG("Reading geometry file \'%s\'.", geometry_file.c_str());

		FileIO::BoostXmlGmlInterface gml_reader(geometries);
		gml_reader.readFile(geometry_file);
	}

	// Process variables
	std::vector<OGS::ProcessVariable> process_variables;

	{
		DBUG("Reading process variables:")
		ConfigTree const& process_variables_config =
			project_config.get_child("process_variables");

		for (auto it : process_variables_config)
		{
			ConfigTree const& var_config = it.second;
			process_variables.emplace_back(var_config, *mesh, geometries);
		}
	}

	// Processes
	std::vector<ProcessLib::Process*> processes;

	{
		ConfigTree processes_config;
		processes_config = project_config.get_child("processes");
		//write_info(std::cout, processes_config);

		DBUG("Reading processes:\n");
		for (auto pc_it : processes_config)
		{
			ConfigTree const& process_config = pc_it.second;

			if (process_config.get<std::string>("type") == "GROUNDWATER_FLOW")
				processes.push_back(new ProcessLib::GroundwaterFlowProcess(
							*mesh, process_variables, process_config));
		}

	}


	std::remove_if(processes.begin(), processes.end(),
			[](ProcessLib::Process* p) { delete p; return true; });
	delete mesh;

	delete fmt;
	delete logog_cout;
	LOGOG_SHUTDOWN();

	return 0;
}
