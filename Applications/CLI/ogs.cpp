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

#ifdef USE_MPI
#include <mpi.h>
#endif

#ifdef USE_PETSC
#include <petsc.h>
#endif

// ThirdParty/tclap
#include "tclap/CmdLine.h"

// BaseLib
#include "BaseLib/BuildInfo.h"
#include "BaseLib/FileTools.h"

#include "Applications/ApplicationsLib/LogogSetup.h"
#include "Applications/ApplicationsLib/ProjectData.h"

#include "ProcessLib/NumericsConfig.h"

void solveProcesses(ProjectData &project)
{
	INFO("Solve processes.");

	std::string const out_pref = project.getOutputFilePrefix();

	auto &time_stepper = project.getTimeStepper();

	while (time_stepper.next())  // skips zeroth timestep, but OK since end of
	                             // first timestep is after first delta t
	{
		const auto dt = time_stepper.getTimeStep().dt();
		const auto current_time = time_stepper.getTimeStep().current();
		const auto timestep = time_stepper.getTimeStep().steps();

		INFO("=================== timestep %i === %g s ===================",
		     timestep, current_time);

		bool accepted = true;

		unsigned i = 0;  // process counter, used to distinguish output files
		for (auto p = project.processesBegin(); p != project.processesEnd();
		     ++p)
		{
			accepted = accepted && (*p)->solve(dt);

			if (!accepted)
			{
				ERR("Timestep has not been accepted. Aborting.");
				break;
			}

			std::string const output_file_name =
			    out_pref + "_pcs_" + std::to_string(i) + "_ts_" +
			    std::to_string(timestep) + ".vtu";
			(*p)->postTimestep(output_file_name, timestep);

			++i;
		}

		if (!accepted)
			break;
	}
}

int main(int argc, char *argv[])
{
	using ConfigTree = boost::property_tree::ptree;

#ifdef USE_MPI
	MPI_Init(&argc, &argv);
#endif

#ifdef USE_PETSC
	char help[] = "ogs6 with PETSc \n";
	PetscInitialize(&argc, &argv, nullptr, help);
#endif

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

	ApplicationsLib::LogogSetup logog_setup;
	
#ifdef USE_LIS
	lis_initialize(&argc, &argv);
#endif

	// Project's configuration
	ConfigTree project_config;

	read_xml(project_arg.getValue(), project_config,
			boost::property_tree::xml_parser::no_comments
			 | boost::property_tree::xml_parser::trim_whitespace);
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

	solveProcesses(project);

#ifdef USE_PETSC
	PetscFinalize();
#endif
#ifdef USE_MPI
	MPI_Finalize();
#endif

#ifdef USE_LIS
	lis_finalize();
#endif

	return 0;
}
