/**
 * \file
 * \author Thomas Fischer
 * \author Wenqing Wang
 * \date   2012-05-10, 2014-10.10
 * \brief  Definition of the RunTime class.
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#ifdef USE_PETSC
#include <mpi.h>
#else
#include <chrono>
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
#ifdef USE_PETSC
        _start_time = MPI_Wtime();
#else
        _start_time = std::chrono::system_clock::now();
#endif
    }

    /// Get the elapsed time in seconds.
    double elapsed() const
    {
#ifdef USE_PETSC
        return MPI_Wtime() - _start_time;
#else
        using namespace std::chrono;
        return duration<double>(system_clock::now() - _start_time).count();
#endif
    }

private:
#ifdef USE_PETSC
    double _start_time = std::numeric_limits<double>::quiet_NaN();
#else
    std::chrono::time_point<std::chrono::system_clock> _start_time;
#endif
};

}  // end namespace BaseLib
