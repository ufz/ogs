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

/// Initialize MPI, PETSc, LIS or any other database from third party
/// packages.
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

/// Release MPI related memory in project, and finalize MPI, PETSc,
/// LIS or any other database from third party
/// packages.
void OgsFinalize(ProjectData &project)
{
#ifdef USE_PETSC
	// Since the global matrix, vector and linear equation in ProjectData
	// are defined as smarter point type (unique_ptr or might be shared_ptr)
	// variables, their memory occupations are automatically released at
	// the end of ProjectData terminated, i.e. the end of the main program.
	// However if  the global matrix, vector and linear equation are created
	// under MPI environment, their memory occupations must be released before
	// calling of MPI_Finalize, PetscFinalize. That is why the following loop is needed.
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

#endif  // OGS_INIT_FINALIZE_H
