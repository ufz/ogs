/**
 * \date   2014-08-04
 * \brief  Implementation of OpenGeoSys simulation application
 *
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */


// ThirdParty/tclap
#include "tclap/CmdLine.h"

// BaseLib
#include "BaseLib/BuildInfo.h"
#include "BaseLib/ConfigTree.h"
#include "BaseLib/ConfigTreeNew.h"
#include "BaseLib/FileTools.h"

#include "Applications/ApplicationsLib/LinearSolverLibrarySetup.h"
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
	// Parse CLI arguments.
	TCLAP::CmdLine cmd("OpenGeoSys-6 software.\n"
			"Copyright (c) 2012-2016, OpenGeoSys Community "
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
	ApplicationsLib::LinearSolverLibrarySetup linear_solver_library_setup(
	    argc, argv);

	// Project's configuration
	BaseLib::ConfigTree project_config =
	    BaseLib::read_xml_config(project_arg.getValue());

	BaseLib::ConfigTreeNew conf(project_config.get_child("OpenGeoSysProject"));
	ProjectData project(conf, BaseLib::extractPath(project_arg.getValue()));

	// Create processes.
	project.buildProcesses<GlobalSetupType>();

	INFO("Initialize processes.");
	for (auto p_it = project.processesBegin(); p_it != project.processesEnd(); ++p_it)
	{
		(*p_it)->initialize();
	}

	std::string const output_file_name(project.getOutputFilePrefix() + ".vtu");

	solveProcesses(project);

	return 0;
}
