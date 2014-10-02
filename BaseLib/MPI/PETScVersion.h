/*!
  \file PETScVersion.h
  \author Wenqing Wang
  \date   2014.09
  \brief  Concatenate PETSc version macros together.

  \copyright
  Copyright (c) 2014, OpenGeoSys Community (http://www.opengeosys.org)
             Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license

*/

#ifndef PETSC_VERSION_H
#define PETSC_VERSION_H

#include <petscversion.h>

#define PETSC_VERSION_NUMBER PETSC_VERSION_MAJOR*1000+PETSC_VERSION_MINOR*10

#endif

