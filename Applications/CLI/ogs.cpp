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
#include "BaseLib/RunTime.h"
#include "InfoLib/CMakeInfo.h"
#include "InfoLib/GitInfo.h"
#include "NumLib/NumericsConfig.h"
#include "ProcessLib/TimeLoop.h"

#ifdef OGS_USE_PYTHON
#include "ogs_embedded_python.h"
#endif

int main(int argc, char* argv[])
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

    TCLAP::ValueArg<std::string> reference_path_arg(
        "r", "reference",
        "Run output result comparison after successful simulation comparing to "
        "all files in the given path. This requires test definitions to be "
        "present in the project file.",
        false, "", "PATH");
    cmd.add(reference_path_arg);

    TCLAP::UnlabeledValueArg<std::string> project_arg(
        "project-file",
        "Path to the ogs6 project file.",
        true,
        "",
        "PROJECT_FILE");
    cmd.add(project_arg);

    TCLAP::MultiArg<std::string> xml_patch_files(
        "p", "xml-patch",
        "the xml patch file(s) which is (are) applied (in the given order) to "
        "the PROJECT_FILE",
        false, "");
    cmd.add(xml_patch_files);

    TCLAP::ValueArg<std::string> outdir_arg("o", "output-directory",
                                            "the output directory to write to",
                                            false, "", "PATH");
    cmd.add(outdir_arg);

    TCLAP::ValueArg<std::string> log_level_arg(
        "l", "log-level",
        "the verbosity of logging messages: none, error, warn, info, debug, "
        "all",
        false,
#ifdef NDEBUG
        "info",
#else
        "all",
#endif
        "LOG_LEVEL");
    cmd.add(log_level_arg);

    TCLAP::SwitchArg nonfatal_arg("",
                                  "config-warnings-nonfatal",
                                  "warnings from parsing the configuration "
                                  "file will not trigger program abortion");
    cmd.add(nonfatal_arg);

    TCLAP::SwitchArg unbuffered_cout_arg("", "unbuffered-std-out",
                                         "use unbuffered standard output");
    cmd.add(unbuffered_cout_arg);

#ifndef _WIN32  // TODO: On windows floating point exceptions are not handled
                // currently
    TCLAP::SwitchArg enable_fpe_arg("", "enable-fpe",
                                    "enables floating point exceptions");
    cmd.add(enable_fpe_arg);
#endif  // _WIN32

    cmd.parse(argc, argv);

    // deactivate buffer for standard output if specified
    if (unbuffered_cout_arg.isSet())
    {
        std::cout.setf(std::ios::unitbuf);
    }

    BaseLib::setConsoleLogLevel(log_level_arg.getValue());
    spdlog::set_pattern("%^%l:%$ %v");
    spdlog::set_error_handler(
        [](const std::string& msg)
        {
            std::cerr << "spdlog error: " << msg << std::endl;
            OGS_FATAL("spdlog logger error occurred.");
        });

    INFO("This is OpenGeoSys-6 version {:s}.",
         GitInfoLib::GitInfo::ogs_version);

#ifndef _WIN32  // On windows this command line option is not present.
    // Enable floating point exceptions
    if (enable_fpe_arg.isSet())
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

    try
    {
        bool solver_succeeded = false;
        {
            ApplicationsLib::LinearSolverLibrarySetup
                linear_solver_library_setup(argc, argv);
#if defined(USE_PETSC)
            vtkSmartPointer<vtkMPIController> controller =
                vtkSmartPointer<vtkMPIController>::New();
            controller->Initialize(&argc, &argv, 1);
            vtkMPIController::SetGlobalController(controller);

            {  // Can be called only after MPI_INIT.
                int mpi_rank;
                MPI_Comm_rank(PETSC_COMM_WORLD, &mpi_rank);
                spdlog::set_pattern(fmt::format("[{}] %^%l:%$ %v", mpi_rank));
            }
#endif

            run_time.start();

            auto project_config = BaseLib::makeConfigTree(
                project_arg.getValue(), !nonfatal_arg.getValue(),
                "OpenGeoSysProject", xml_patch_files.getValue());

            BaseLib::setProjectDirectory(
                BaseLib::extractPath(project_arg.getValue()));

            ProjectData project(*project_config,
                                BaseLib::getProjectDirectory(),
                                outdir_arg.getValue());

            if (!reference_path_arg.isSet())
            {  // Ignore the test_definition section.
                project_config->ignoreConfigParameter("test_definition");
            }
            else
            {
                test_definition =
                    std::make_unique<ApplicationsLib::TestDefinition>(
                        //! \ogs_file_param{prj__test_definition}
                        project_config->getConfigSubtree("test_definition"),
                        reference_path_arg.getValue(),
                        outdir_arg.getValue());
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
            project_config->ignoreConfigParameter("insitu");
#endif

            INFO("Initialize processes.");
            for (auto& p : project.getProcesses())
            {
                p->initialize();
            }

            // Check intermediately that config parsing went fine.
            project_config.checkAndInvalidate();
            BaseLib::ConfigTree::assertNoSwallowedErrors();

            INFO("Solve processes.");

            auto& time_loop = project.getTimeLoop();
            time_loop.initialize();
            solver_succeeded = time_loop.loop();

#ifdef USE_INSITU
            if (isInsituConfigured)
                InSituLib::Finalize();
#endif
            INFO("[time] Execution took {:g} s.", run_time.elapsed());

#if defined(USE_PETSC)
            controller->Finalize(1);
#endif
        }  // This nested scope ensures that everything that could possibly
           // possess a ConfigTree is destructed before the final check below is
           // done.

        BaseLib::ConfigTree::assertNoSwallowedErrors();

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
