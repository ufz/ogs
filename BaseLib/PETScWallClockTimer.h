/*!
  \file PETScWallClockTimer.h
  \author Wenqing Wang
  \date   2014.08
  \brief  Declare a class to record wall clock time in computation with PETSc.

  \copyright
  Copyright (c) 2014, OpenGeoSys Community (http://www.opengeosys.org)
             Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license

*/

#ifndef PETSC_WALL_CLOCK_TIMER_H
#define PETSC_WALL_CLOCK_TIMER_H

#include <petsctime.h>

namespace BaseLib
{

/// Record wall clock time for computations with PETSc.
class PETScWallClockTimer
{
    public:
        /// Record the start time.
        void start()
        {
#if (PETSC_VERSION_NUMBER > 3020)
            PetscTime(&_start_time);
#else
            PetscGetTime(&_start_time);
#endif
        }

        /// Return the elapsed time when this function is called.
        PetscLogDouble elapsed() const
        {
            PetscLogDouble current_time;
#if (PETSC_VERSION_NUMBER > 3020)
            PetscTime(&current_time);
#else
            PetscGetTime(&current_time);
#endif
            return current_time - _start_time;
        }

    private:
        /// Start time.
        PetscLogDouble _start_time = 0.;
};

} // end namespace BaseLib

#endif

