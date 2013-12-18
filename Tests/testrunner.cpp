/**
 * \file
 * \author Lars Bilke
 * \date   2012-04-29
 * \brief  GTest test executables main function.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

// ** INCLUDES **
#include "Configure.h"
#include "gtest/gtest.h"
#include "logog/include/logog.hpp"

#ifdef USE_LIS
#include "lis.h"
#endif

#ifdef OGS_USE_PETSC
#include "petscksp.h"
#endif

#ifdef OGS_USE_BOOSTMPI // boost must be included after petsc headers
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#else
#include "MathLib/LinAlg/PETSc/InfoMPI.h"
#endif

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

#if defined(OGS_USE_PETSC)
        char help[] = "ogs6 with PETSc \n";
        PetscInitialize(&argc,&argv,(char *)0,help);

        int mrank, msize;
        MPI_Comm_rank(PETSC_COMM_WORLD, &mrank);
        MPI_Comm_size(PETSC_COMM_WORLD, &msize);

#ifdef OGS_USE_BOOSTMPI // demo
        boost::mpi::communicator petsc( PETSC_COMM_WORLD, boost::mpi::comm_attach);
#else
        BaseLib::InfoMPI::setSizeRank(msize, mrank);
#endif
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

#if defined(OGS_USE_PETSC)
        PetscFinalize();
#endif


    } // make sure no logog objects exist when LOGOG_SHUTDOWN() is called.
    LOGOG_SHUTDOWN();

    return ret;
}
