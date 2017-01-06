/**
 * \file
 * \author Lars Bilke
 * \date   2012-04-29
 * \brief  GTest test executables main function.
 *
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

// ** INCLUDES **
#include <clocale>

#include "gtest/gtest.h"

#include "Applications/ApplicationsLib/LogogSetup.h"
#include "Applications/ApplicationsLib/LinearSolverLibrarySetup.h"
#include "NumLib/DOF/GlobalMatrixProviders.h"
#include "BaseLib/BuildInfo.h"
#include "BaseLib/TemplateLogogFormatterSuppressedGCC.h"

#ifdef OGS_BUILD_GUI
#include <QCoreApplication>
#endif

/// Implementation of the googletest testrunner
int main(int argc, char* argv[])
{
    std::string logLevel("all");
    for (int i = 1; i < argc; i++)
    {
        if(i + 1 == argc)
            break;
        if(std::strcmp(argv[i], "-l") == 0)
            logLevel = argv[i + 1];
    }

    setlocale(LC_ALL, "C");
#ifdef OGS_BUILD_GUI
    QCoreApplication app(argc, argv, false);
#endif

    // Attention: Order matters!
    // logog_setup must be created first, then the linear_solver_library_setup,
    // because destruction order is the reverse of the creation order and the
    // destructor of linear_solver_library_setup might print log messages.
    // The methods on logog_setup must be called after the construction of
    // linear_solver_library_setup since they require, e.g., that MPI_Init()
    // has been called before.

    ApplicationsLib::LogogSetup logog_setup;

    ApplicationsLib::LinearSolverLibrarySetup linear_solver_library_setup(
                argc, argv);

    logog_setup.setFormatter(std::unique_ptr<BaseLib::TemplateLogogFormatterSuppressedGCC
        <TOPIC_LEVEL_FLAG | TOPIC_FILE_NAME_FLAG | TOPIC_LINE_NUMBER_FLAG> >
            (new BaseLib::TemplateLogogFormatterSuppressedGCC
            <TOPIC_LEVEL_FLAG | TOPIC_FILE_NAME_FLAG | TOPIC_LINE_NUMBER_FLAG>()));
    logog_setup.setLevel(logLevel);

    try
    {
        // start google test
        testing::InitGoogleTest ( &argc, argv );
        return RUN_ALL_TESTS();
    }
    catch (char* e)
    {
        ERR(e);
        return 1;
    }
    catch (std::exception& e)
    {
        ERR(e.what());
        return 1;
    }
    catch (...)
    {
        ERR("Unknown exception occurred!");
        return 1;
    }

    return 0;
}
