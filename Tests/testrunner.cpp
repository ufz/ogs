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

#ifdef USE_MPI
#include <mpi.h>
#endif

#ifdef USE_LIS
#include <lis.h>
#endif

#ifdef USE_PETSC
#include <petscksp.h>
#endif

#include "Applications/ApplicationsLib/LogogSetup.h"
#include "BaseLib/TemplateLogogFormatterSuppressedGCC.h"
#ifdef OGS_BUILD_GUI
#include <QApplication>
#endif

/// Implementation of the googletest testrunner
int main(int argc, char* argv[])
{
    setlocale(LC_ALL, "C");
#ifdef OGS_BUILD_GUI
    QApplication app(argc, argv, false);
#endif
    int ret = 0;
#ifdef USE_MPI
    MPI_Init(&argc, &argv);
#endif
    ApplicationsLib::LogogSetup logog_setup;
    logog_setup.SetFormatter(std::unique_ptr<BaseLib::TemplateLogogFormatterSuppressedGCC
        <TOPIC_LEVEL_FLAG | TOPIC_FILE_NAME_FLAG | TOPIC_LINE_NUMBER_FLAG> >
            (new BaseLib::TemplateLogogFormatterSuppressedGCC
            <TOPIC_LEVEL_FLAG | TOPIC_FILE_NAME_FLAG | TOPIC_LINE_NUMBER_FLAG>()));


#ifdef USE_PETSC
    char help[] = "ogs6 with PETSc \n";
    PetscInitialize(&argc,&argv,(char *)0,help);
#endif

    try
    {
        // initialize libraries which will be used while testing
#ifdef USE_LIS
        lis_initialize(&argc, &argv);
#endif
        // start google test
        testing::InitGoogleTest ( &argc, argv );
        ret = RUN_ALL_TESTS();
    }
    catch (char* e)
    {
        ERR(e);
    }
    catch (std::exception& e)
    {
        ERR(e.what());
    }
    catch (...)
    {
        ERR("Unknown exception occurred!");
    }
    // finalize libraries
#ifdef USE_LIS
    lis_finalize();
#endif

#ifdef USE_PETSC
    PetscFinalize();
#endif

#ifdef USE_MPI
    MPI_Finalize();
#endif

    return ret;
}
