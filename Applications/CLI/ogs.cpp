/**
 * \date   2014-08-04
 * \brief  Implementation of OpenGeoSys simulation application
 *
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <tclap/CmdLine.h>
#include <chrono>

#ifndef _WIN32
#ifdef __APPLE__
#include <xmmintrin.h>
#else
#include <cfenv>
#endif  // __APPLE__
#endif  // _WIN32

#ifdef USE_PETSC
#include <vtkMPIController.h>
#include <vtkSmartPointer.h>
#endif

// BaseLib
#include "BaseLib/BuildInfo.h"
#include "BaseLib/ConfigTreeUtil.h"
#include "BaseLib/DateTools.h"
#include "BaseLib/FileTools.h"
#include "BaseLib/RunTime.h"
#include "BaseLib/TemplateLogogFormatterSuppressedGCC.h"

#include "Applications/ApplicationsLib/LinearSolverLibrarySetup.h"
#include "Applications/ApplicationsLib/LogogSetup.h"
#include "Applications/ApplicationsLib/ProjectData.h"
#include "Applications/ApplicationsLib/TestDefinition.h"
#include "Applications/InSituLib/Adaptor.h"
#include "ProcessLib/TimeLoop.h"

#include "NumLib/NumericsConfig.h"

#ifdef OGS_USE_PYTHON
#include "ogs_embedded_python.h"
#endif

int main(int argc, char* argv[])
{
    // Parse CLI arguments.
    TCLAP::CmdLine cmd(
        "OpenGeoSys-6 software.\n"
        "Copyright (c) 2012-2019, OpenGeoSys Community "
        "(http://www.opengeosys.org) "
        "Distributed under a Modified BSD License. "
        "See accompanying file LICENSE.txt or "
        "http://www.opengeosys.org/project/license\n"
        "version: " +
            BaseLib::BuildInfo::ogs_version + "\n" +
        "CMake arguments: " +
            BaseLib::BuildInfo::cmake_args,
        ' ',
        BaseLib::BuildInfo::ogs_version);

    TCLAP::ValueArg<std::string> reference_path_arg(
        "r", "reference",
        "Run output result comparison after successful simulation comparing to "
        "all files in the given path.",
        false, "", "PATH");
    cmd.add(reference_path_arg);

    TCLAP::UnlabeledValueArg<std::string> project_arg(
        "project-file",
        "Path to the ogs6 project file.",
        true,
        "",
        "PROJECT_FILE");
    cmd.add(project_arg);

    TCLAP::ValueArg<std::string> outdir_arg("o", "output-directory",
                                            "the output directory to write to",
                                            false, "", "PATH");
    cmd.add(outdir_arg);

    TCLAP::ValueArg<std::string> log_level_arg("l", "log-level",
                                               "the verbosity of logging "
                                               "messages: none, error, warn, "
                                               "info, debug, all",
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

    ApplicationsLib::LogogSetup logog_setup;
    logog_setup.setLevel(log_level_arg.getValue());

    INFO("This is OpenGeoSys-6 version %s.",
         BaseLib::BuildInfo::ogs_version.c_str());

#ifndef _WIN32  // On windows this command line option is not present.
    // Enable floating point exceptions
    if (enable_fpe_arg.isSet())
    {
#ifdef __APPLE__
        _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_INVALID);
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
        INFO("OGS started on %s.", time_str.c_str());
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

            logog_setup.setFormatter(
                std::make_unique<BaseLib::TemplateLogogFormatterSuppressedGCC<
                    TOPIC_LEVEL_FLAG | TOPIC_FILE_NAME_FLAG |
                    TOPIC_LINE_NUMBER_FLAG>>());
#endif
            run_time.start();

            auto project_config = BaseLib::makeConfigTree(
                project_arg.getValue(), !nonfatal_arg.getValue(),
                "OpenGeoSysProject");

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

            BaseLib::ConfigTree::assertNoSwallowedErrors();

            BaseLib::ConfigTree::assertNoSwallowedErrors();

            INFO("Solve processes.");

            auto& time_loop = project.getTimeLoop();
            solver_succeeded = time_loop.loop();

#ifdef USE_INSITU
            if (isInsituConfigured)
                InSituLib::Finalize();
#endif
            INFO("[time] Execution took %g s.", run_time.elapsed());

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
        INFO("OGS terminated on %s.", time_str.c_str());
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
