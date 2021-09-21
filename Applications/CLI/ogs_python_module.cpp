/**
 * \brief  Implementation of OpenGeoSys simulation application python module
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <spdlog/spdlog.h>
#include <tclap/CmdLine.h>

#include "CommandLineArgumentParser.h"

#include "Applications/ApplicationsLib/Simulation.h"

#include "Applications/ApplicationsLib/TestDefinition.h"
#include "BaseLib/DateTools.h"
#include "BaseLib/Error.h"
#include "BaseLib/Logging.h"
#include "BaseLib/RunTime.h"
#include "InfoLib/GitInfo.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#ifdef OGS_USE_PYTHON
#include "ogs_embedded_python.h"
#endif

void initOGS(std::vector<std::string>& argv_str)
{
    std::cout << "begin C++ world " << std::endl;

    int argc = argv_str.size();
    char **argv = new char*[argc];
    for (int i = 0; i < argc; ++i)
    {
        argv[i] = argv_str[i].data();
    }

    CommandLineArgumentParser cli_args(argc, argv);
    BaseLib::initOGSLogger(cli_args.log_level);

    INFO("This is OpenGeoSys-6 version {:s}.",
         GitInfoLib::GitInfo::ogs_version);

    BaseLib::RunTime run_time;

    {
        auto const start_time = std::chrono::system_clock::now();
        auto const time_str = BaseLib::formatDate(start_time);
        INFO("OGS started on {:s}.", time_str);
    }

    std::optional<ApplicationsLib::TestDefinition> test_definition{
        std::nullopt};
    auto ogs_status = EXIT_SUCCESS;

    try
    {
        Simulation simulation(argc, argv);
        run_time.start();
        simulation.initializeDataStructures(
            std::move(cli_args.project), std::move(cli_args.xml_patch_file_names),
            cli_args.reference_path_is_set, std::move(cli_args.reference_path),
            cli_args.nonfatal, std::move(cli_args.outdir),
            std::move(cli_args.mesh_dir), cli_args.write_prj);
        bool solver_succeeded = simulation.executeSimulation();
        test_definition = simulation.getTestDefinition();

        INFO("[time] Execution took {:g} s.", run_time.elapsed());
        ogs_status = solver_succeeded ? EXIT_SUCCESS : EXIT_FAILURE;
    }
    catch (std::exception& e)
    {
        ERR(e.what());
        ogs_status = EXIT_FAILURE;
    }

    {
        auto const end_time = std::chrono::system_clock::now();
        auto const time_str = BaseLib::formatDate(end_time);
        INFO("OGS terminated on {:s}.", time_str);
    }

    std::cout << "OpenGeoSys module initOGS was called" << std::endl;
    std::copy(argv_str.begin(), argv_str.end(),
              std::ostream_iterator<std::string>(std::cout, " "));
    std::cout << std::endl;
    std::cout << "end C++ world " << std::endl;
}

/// python module name is OpenGeoSys
PYBIND11_MODULE(OpenGeoSys, m)
{
    m.doc() = "pybind11 ogs example plugin";
    m.def("initOGS", &initOGS, "init OGS");
}
