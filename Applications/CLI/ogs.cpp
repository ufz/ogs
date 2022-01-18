/**
 * \date   2014-08-04
 * \brief  Implementation of OpenGeoSys simulation application
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

#include <chrono>
#include <sstream>

#ifndef _WIN32
#ifdef __APPLE__
#ifdef __SSE__
#include <xmmintrin.h>
#endif  // __SSE__
#else
#include <cfenv>
#endif  // __APPLE__
#endif  // _WIN32

#ifdef USE_PETSC
#include <vtkMPIController.h>
#include <vtkSmartPointer.h>
#endif

// BaseLib
#include "Applications/ApplicationsLib/LinearSolverLibrarySetup.h"
#include "Applications/ApplicationsLib/ProjectData.h"
#include "Applications/ApplicationsLib/TestDefinition.h"
#include "Applications/InSituLib/Adaptor.h"
#include "BaseLib/ConfigTreeUtil.h"
#include "BaseLib/DateTools.h"
#include "BaseLib/FileTools.h"
#include "BaseLib/PrjProcessing.h"
#include "BaseLib/RunTime.h"
#include "InfoLib/CMakeInfo.h"
#include "InfoLib/GitInfo.h"
#include "NumLib/NumericsConfig.h"
#include "ProcessLib/TimeLoop.h"

#ifdef OGS_USE_PYTHON
#include "ogs_embedded_python.h"
#endif

#ifndef _WIN32  // On windows this command line option is not present.
void enableFloatingPointExceptions()
    {
#ifdef __APPLE__
#ifdef __SSE__
        _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_INVALID);
#endif  // __SSE__
#else
        feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#endif  // __APPLE__
    }
#endif  // _WIN32
void setConsoleLogLevel(TCLAP::ValueArg<std::string> const& log_level_arg)
{
    BaseLib::setConsoleLogLevel(log_level_arg.getValue());
    spdlog::set_pattern("%^%l:%$ %v");
    spdlog::set_error_handler(
        [](const std::string& msg)
        {
            std::cerr << "spdlog error: " << msg << std::endl;
            OGS_FATAL("spdlog logger error occurred.");
        });
}

struct CommandLineArguments final
{
    CommandLineArguments(int argc, char* argv[])
    {
        // Parse CLI arguments.
        TCLAP::CmdLine cmd(
            "OpenGeoSys-6 software.\n"
            "Copyright (c) 2012-2022, OpenGeoSys Community "
            "(http://www.opengeosys.org) "
            "Distributed under a Modified BSD License. "
            "See accompanying file LICENSE.txt or "
            "http://www.opengeosys.org/project/license\n"
            "version: " +
                GitInfoLib::GitInfo::ogs_version + "\n" +
                "CMake arguments: " + CMakeInfoLib::CMakeInfo::cmake_args,
            ' ',
            GitInfoLib::GitInfo::ogs_version + "\n\n" +
                "CMake arguments: " + CMakeInfoLib::CMakeInfo::cmake_args);

        TCLAP::ValueArg<std::string> log_level_arg(
            "l", "log-level",
            "the verbosity of logging messages: none, error, warn, info, "
            "debug, "
            "all",
            false,
#ifdef NDEBUG
            "info",
#else
            "all",
#endif
            "LOG_LEVEL");

#ifndef _WIN32  // TODO: On windows floating point exceptions are not handled
                // currently
        TCLAP::SwitchArg enable_fpe_arg("", "enable-fpe",
                                        "enables floating point exceptions");
#endif  // _WIN32
        TCLAP::SwitchArg unbuffered_cout_arg("", "unbuffered-std-out",
                                             "use unbuffered standard output");

        TCLAP::ValueArg<std::string> reference_path_arg(
            "r", "reference",
            "Run output result comparison after successful simulation "
            "comparing to all files in the given path. This requires test "
            "definitions to be present in the project file.",
            false, "", "PATH");

        TCLAP::UnlabeledValueArg<std::string> project_arg(
            "project-file",
            "Path to the ogs6 project file.",
            true,
            "",
            "PROJECT_FILE");

        TCLAP::MultiArg<std::string> xml_patch_files_arg(
            "p", "xml-patch",
            "the xml patch file(s) which is (are) applied (in the given order) "
            "to the PROJECT_FILE",
            false, "");

        TCLAP::ValueArg<std::string> outdir_arg(
            "o", "output-directory", "the output directory to write to", false,
            "", "PATH");

        TCLAP::SwitchArg nonfatal_arg("",
                                      "config-warnings-nonfatal",
                                      "warnings from parsing the configuration "
                                      "file will not trigger program abortion");
        cmd.add(reference_path_arg);
        cmd.add(project_arg);
        cmd.add(xml_patch_files_arg);
        cmd.add(outdir_arg);
        cmd.add(log_level_arg);
        cmd.add(nonfatal_arg);
        cmd.add(unbuffered_cout_arg);
#ifndef _WIN32  // TODO: On windows floating point exceptions are not handled
                // currently
        cmd.add(enable_fpe_arg);
#endif  // _WIN32

        cmd.parse(argc, argv);

        reference_path = reference_path_arg.getValue();
        reference_path_is_set = reference_path_arg.isSet();
        project = project_arg.getValue();
        xml_patch_file_names = xml_patch_files_arg.getValue();
        outdir = outdir_arg.getValue();
        nonfatal = nonfatal_arg.getValue();

        // deactivate buffer for standard output if specified
        if (unbuffered_cout_arg.isSet())
        {
            std::cout.setf(std::ios::unitbuf);
        }

        setConsoleLogLevel(log_level_arg);
#ifndef _WIN32  // TODO: On windows floating point exceptions are not handled
        if (enable_fpe_arg.isSet())
        {
            enableFloatingPointExceptions();
        }
#endif  // _WIN32
    }

    std::string reference_path;
    std::string project;
    std::vector<std::string> xml_patch_file_names;
    std::string outdir;
    bool nonfatal;
    bool reference_path_is_set;
};

int main(int argc, char* argv[])
{
    CommandLineArguments cli_arg(argc, argv);

    INFO("This is OpenGeoSys-6 version {:s}.",
         GitInfoLib::GitInfo::ogs_version);

#ifdef OGS_USE_PYTHON
    pybind11::scoped_interpreter guard = ApplicationsLib::setupEmbeddedPython();
    (void)guard;
#endif

    BaseLib::RunTime run_time;

    {
        auto const start_time = std::chrono::system_clock::now();
        auto const time_str = BaseLib::formatDate(start_time);
        INFO("OGS started on {:s}.", time_str);
    }

    std::unique_ptr<ApplicationsLib::TestDefinition> test_definition;
    auto ogs_status = EXIT_SUCCESS;

#if defined(USE_PETSC)
    vtkSmartPointer<vtkMPIController> controller;
#endif
    std::unique_ptr<ProjectData> project;
    try
    {
        ApplicationsLib::LinearSolverLibrarySetup linear_solver_library_setup(
            argc, argv);
#if defined(USE_PETSC)
        controller = vtkSmartPointer<vtkMPIController>::New();
        controller->Initialize(&argc, &argv, 1);
        vtkMPIController::SetGlobalController(controller);

        {  // Can be called only after MPI_INIT.
            int mpi_rank;
            MPI_Comm_rank(PETSC_COMM_WORLD, &mpi_rank);
            spdlog::set_pattern(fmt::format("[{}] %^%l:%$ %v", mpi_rank));
        }
#endif

        run_time.start();

        {
            auto project_config = BaseLib::makeConfigTree(
                cli_arg.project, !cli_arg.nonfatal, "OpenGeoSysProject",
                cli_arg.xml_patch_file_names);

            BaseLib::setProjectDirectory(BaseLib::extractPath(cli_arg.project));

            if (!cli_arg.reference_path_is_set)
            {  // Ignore the test_definition section.
                project_config.ignoreConfigParameter("test_definition");
            }
            else
            {
                test_definition =
                    std::make_unique<ApplicationsLib::TestDefinition>(
                        //! \ogs_file_param{prj__test_definition}
                        project_config.getConfigSubtree("test_definition"),
                        cli_arg.reference_path, cli_arg.outdir);
                if (test_definition->numberOfTests() == 0)
                {
                    OGS_FATAL(
                        "No tests were constructed from the test definitions, "
                        "but reference solutions path was given.");
                }

                INFO("Cleanup possible output files before running ogs.");
                BaseLib::removeFiles(test_definition->getOutputFiles());
            }
#ifdef USE_INSITU
            auto isInsituConfigured = false;
            //! \ogs_file_param{prj__insitu}
            if (auto t = project_config->getConfigSubtreeOptional("insitu"))
            {
                InSituLib::Initialize(
                    //! \ogs_file_param{prj__insitu__scripts}
                    t->getConfigSubtree("scripts"),
                    BaseLib::extractPath(project_arg.getValue()));
                isInsituConfigured = true;
            }
#else
            project_config.ignoreConfigParameter("insitu");
#endif

            project =
                std::make_unique<ProjectData>(project_config,
                                              BaseLib::getProjectDirectory(),
                                              cli_arg.outdir);

            INFO("Initialize processes.");
            for (auto& p : project->getProcesses())
            {
                p->initialize();
            }

            // Check intermediately that config parsing went fine.
            checkAndInvalidate(project_config);
            BaseLib::ConfigTree::assertNoSwallowedErrors();

            auto& time_loop = project->getTimeLoop();
            time_loop.initialize();
        }

        bool solver_succeeded = false;
        {
            INFO("Solve processes.");
            auto& time_loop = project->getTimeLoop();
            solver_succeeded = time_loop.loop();
        }

        {
#ifdef USE_INSITU
            if (isInsituConfigured)
                InSituLib::Finalize();
#endif
            INFO("[time] Execution took {:g} s.", run_time.elapsed());

#if defined(USE_PETSC)
            controller->Finalize(1);
#endif
            // This nested scope ensures that everything that could possibly
            // possess a ConfigTree is destructed before the final check below
            // is done.

            INFO("cleanup:  ProjectData ...");
            auto* project_pointer = project.release();
            delete project_pointer;
            INFO("cleanup:  ProjectData done.");

            //INFO("cleanup:  BaseLib::ConfigTree::assertNoSwallowedErrors()...");
            //BaseLib::ConfigTree::assertNoSwallowedErrors();

            ogs_status = solver_succeeded ? EXIT_SUCCESS : EXIT_FAILURE;
        }
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

    if (ogs_status == EXIT_FAILURE)
    {
        ERR("OGS terminated with error.");
        return EXIT_FAILURE;
    }

    if (test_definition == nullptr)
    {
        // There are no tests, so just exit;
        return ogs_status;
    }

    INFO("");
    INFO("##########################################");
    INFO("# Running tests                          #");
    INFO("##########################################");
    INFO("");
    if (!test_definition->runTests())
    {
        ERR("One of the tests failed.");
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
