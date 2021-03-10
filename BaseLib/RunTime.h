/**
 * \file
 * \author Thomas Fischer
 * \author Wenqing Wang
 * \date   2012-05-10, 2014-10.10
 * \brief  Definition of the RunTime class.
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
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
        start_time_ = MPI_Wtime();
#else
        start_time_ = std::chrono::system_clock::now();
#endif
    }

    /// Get the elapsed time in seconds.
    double elapsed() const
    {
#ifdef USE_PETSC
        return MPI_Wtime() - start_time_;
#else
        using namespace std::chrono;
        return duration<double>(system_clock::now() - start_time_).count();
#endif
    }

private:
#ifdef USE_PETSC
    double start_time_ = std::numeric_limits<double>::quiet_NaN();
#else
    std::chrono::time_point<std::chrono::system_clock> start_time_;
#endif
};

}  // end namespace BaseLib
