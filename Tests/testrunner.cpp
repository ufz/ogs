/**
 * \file
 * \author Lars Bilke
 * \date   2012-04-29
 * \brief  GTest test executables main function.
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

// ** INCLUDES **
#ifdef USE_LIS
#include <lis.h>
#endif

#ifdef USE_PETSC
#include <petscksp.h>
#endif

#include "gtest/gtest.h"
#include "BaseLib/LogogWrapper.h"
//#include "logog/include/logog.hpp"

#include "BaseLib/TemplateLogogFormatterSuppressedGCC.h"
#ifdef QT4_FOUND
#include <QApplication>
#endif

/// Implementation of the googletest testrunner
int main(int argc, char* argv[])
{
#ifdef QT4_FOUND
	QApplication app(argc, argv, false);
#endif
    int ret = 0;
    LOGOG_INITIALIZE();
    {
        logog::Cout out;
        BaseLib::TemplateLogogFormatterSuppressedGCC<TOPIC_LEVEL_FLAG | TOPIC_FILE_NAME_FLAG | TOPIC_LINE_NUMBER_FLAG> custom_format;
        out.SetFormatter(custom_format);

#ifdef USE_PETSC
        char help[] = "ogs6 with PETSc \n";
        PetscInitialize(&argc, &argv, nullptr, help);
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
            BaseLib::err(e);
        }
        catch (std::exception& e)
        {
            BaseLib::err(e.what());
        }
        catch (...)
        {
            BaseLib::err("Unknown exception occurred!");
        }
        // finalize libraries
#ifdef USE_LIS
        lis_finalize();
#endif

#ifdef USE_PETSC
        PetscFinalize();
#endif

    } // make sure no logog objects exist when LOGOG_SHUTDOWN() is called.
    LOGOG_SHUTDOWN();

    return ret;
}
