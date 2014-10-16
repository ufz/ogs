/**
 * \file
 * \author Thomas Fischer
 * \author Wenqing Wang
 * \date   2012-05-10, 2014-10.10
 * \brief  Definition of the RunTime class.
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef RUNTIME_H
#define RUNTIME_H

#if defined(USE_MPI)
#include <mpi.h>
#else
#ifndef _MSC_VER
#include <sys/time.h>
#else
#include <windows.h>
#endif
#endif

namespace BaseLib
{

/// Count the running time.
class RunTime
{
    public:
        /// Start the timer.
        void start()
        {
#if defined(USE_MPI)
            _timer = MPI_Wtime();
#else
#ifndef _MSC_VER
            timeval t;
            gettimeofday(&t, 0);
            _timer = -t.tv_sec - t.tv_usec/1000000.0;
#else
            _timer = timeGetTime();
#endif
#endif
        }

        /// Get the elapsed time after started.
        double elapsed()
        {
#if defined(USE_MPI)
            return _timer - MPI_Wtime();
#else
#ifndef _MSC_VER
            timeval t;
            gettimeofday(&t, 0);
            _timer += t.tv_sec + t.tv_usec/1000000.0;
            return _timer;
#else
            return (timeGetTime() - _timer)/1000.0;
#endif
#endif
        }

    private:
        double _timer = 0.;
};

} // end namespace BaseLib

#endif
