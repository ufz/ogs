/**
 * \file
 * \author Thomas Fischer
 * \author Wenqing Wang
 * \date   2012-05-10, 2014-10.10
 * \brief  Definition of the RunTime class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
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
#ifndef WIN32
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
            _start_time = MPI_Wtime();
#else
#ifndef _MSC_VER
            timeval t;
            gettimeofday(&t, 0);
            _start_time = t.tv_sec + t.tv_usec/1000000.0;
#else
            _start_time = timeGetTime();
#endif
#endif
        }

        /// Get the elapsed time after started.
        double elapsed() const
        {
#if defined(USE_MPI)
            return MPI_Wtime() - _start_time;
#else
#ifndef _MSC_VER
            timeval t;
            gettimeofday(&t, 0);
            return t.tv_sec + t.tv_usec/1000000.0 - _start_time;
#else
            return (timeGetTime() - _start_time)/1000.0;
#endif
#endif
        }

    private:
        double _start_time = 0.;
};

} // end namespace BaseLib

#endif
