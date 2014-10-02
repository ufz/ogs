/*!
  \file WallClockTimer.h
  \author Wenqing Wang
  \date   2014.08
  \brief  Declare a class to record wall clock time in computation with MPI.

  \copyright
  Copyright (c) 2014, OpenGeoSys Community (http://www.opengeosys.org)
             Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license

*/

#ifndef WALL_CLOCK_TIMER_H
#define WALL_CLOCK_TIMER_H

#include <mpi.h>

namespace BaseLib
{

/// Record wall clock time in the computations with MPI
class WallClockTimer
{
    public:
        /// Record the start time
        void start()
        {
            _timer = -MPI_Wtime();
        }

        /// Return the elapsed time when this function is called.
        double elapsed()
        {
            return _timer + MPI_Wtime();
        }

    private:
        double _timer = 0.;
};

} // end namespace BaseLib

#endif

