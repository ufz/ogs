/*!
  \file OgsInitFinalize
  \brief Define two functions for external algebraic packages:
         one to initialize data base or together with MPI,
         another to close the program with the external algebraic packages 
  \copyright
  Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
             Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license
 
 */

#ifndef OGS_INIT_FINALIZE_H
#define OGS_INIT_FINALIZE_H

#ifdef USE_MPI
#include <mpi.h>
#endif

#ifdef USE_PETSC
#include <petsc.h>
#endif

#ifdef USE_LIS
#include <lis.h>
#endif

#include "Applications/ApplicationsLib/ProjectData.h"

namespace BaseLib
{

void OgsInitialize(int argc, char *argv[])
{
#ifdef USE_MPI
	MPI_Init(&argc, &argv);
#endif

#ifdef USE_PETSC
	char help[] = "ogs6 with PETSc \n";
	PetscInitialize(&argc, &argv, nullptr, help);
#endif

#ifdef USE_LIS
	lis_initialize(&argc, &argv);
#endif

	 // Avoid compilation warning in case none of above cases exsiting. 
	(void)argc;
	(void)argv;    
}

void OgsFinalize(ProjectData &project)
{
#ifdef USE_PETSC
	for (auto p_it = project.processesBegin(); p_it != project.processesEnd(); ++p_it)
	{
		(*p_it)->releaseEquationMemory(); // possibly also needed under USE_MPI.
	}
	PetscFinalize();
#endif

#ifdef USE_MPI
	MPI_Finalize();
#endif

#ifdef USE_LIS
	lis_finalize();
#endif

	 // Avoid compilation warning in case that project is not used. 
	(void) project;
}

} // end namespace BaseLib

#endif  // OGS_INIT_FINALIZE_H
