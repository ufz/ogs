/*!
  \file LogogWrapper.h
  \author Wenqing Wang
  \date   2014.10
  \brief  Wrapper for logog macros to enable them for less verbose output
          under MPI

  \copyright
  Copyright (c) 2014, OpenGeoSys Community (http://www.opengeosys.org)
             Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license

*/

#ifndef LOGOG_WRAPPER_H
#define LOGOG_WRAPPER_H

#ifdef USE_MPI
#include <mpi.h>
#endif

#include "logog/include/logog.hpp"

namespace BaseLib
{

/// Wrapped INFO for common usage.
template<typename... T_LOGOG_ARGS> void info(T_LOGOG_ARGS const& ... args)
{
#ifdef USE_MPI
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if(rank == 0)
#endif
        INFO(args...);
};

/// Wrapped WARN for common usage.
template<typename... T_LOGOG_ARGS> void warn(T_LOGOG_ARGS const& ... args)
{
#ifdef USE_MPI
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if(rank == 0)
#endif
        WARN(args...);
};

/// Wrapped WARN1 for common usage.
template<typename... T_LOGOG_ARGS> void warn1(T_LOGOG_ARGS const& ... args)
{
#ifdef USE_MPI
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if(rank == 0)
#endif
        WARN1(args...);
};

/// Wrapped WARN2 for common usage.
template<typename... T_LOGOG_ARGS> void warn2(T_LOGOG_ARGS const& ... args)
{
#ifdef USE_MPI
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if(rank == 0)
#endif
        WARN2(args...);
};

/// Wrapped WARN3 for common usage.
template<typename... T_LOGOG_ARGS> void warn3(T_LOGOG_ARGS const& ... args)
{
#ifdef USE_MPI
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if(rank == 0)
#endif
        WARN3(args...);
};

/// Wrapped DBUG for common usage.
template<typename... T_LOGOG_ARGS> void dbug(T_LOGOG_ARGS const& ... args)
{
#ifdef USE_MPI
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if(rank == 0)
#endif
        DBUG(args...);
};

/// Wrapped ERR for common usage.
template<typename... T_LOGOG_ARGS> void err(T_LOGOG_ARGS const& ... args)
{
#ifdef USE_MPI
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if(rank == 0)
#endif
        ERR(args...);
};

/// Wrapped ALERT for common usage.
template<typename... T_LOGOG_ARGS> void alert(T_LOGOG_ARGS const& ... args)
{
#ifdef USE_MPI
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if(rank == 0)
#endif
        ALERT(args...);
};

/// Wrapped CRITICAL for common usage.
template<typename... T_LOGOG_ARGS> void critical(T_LOGOG_ARGS const& ... args)
{
#ifdef USE_MPI
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if(rank == 0)
#endif
        CRITICAL(args...);
}

/// Wrapped EMERGENCY for common usage.
template<typename... T_LOGOG_ARGS> void emergency(T_LOGOG_ARGS const& ... args)
{
#ifdef USE_MPI
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if(rank == 0)
#endif
        EMERGENCY(args...);
};

#ifdef USE_MPI
/// Wrapped INFO for specified commnicator.
template<typename... T_LOGOG_ARGS> void info(MPI_Comm comm, T_LOGOG_ARGS const& ... args)
{
    int rank;
    MPI_Comm_rank(comm, &rank);
    if(rank == 0)
        INFO(args...);
};

/// Wrapped WARN for specified commnicator.
template<typename... T_LOGOG_ARGS> void warn(MPI_Comm comm, T_LOGOG_ARGS const& ... args)
{
    int rank;
    MPI_Comm_rank(comm, &rank);
    if(rank == 0)
        WARN(args...);
};

/// Wrapped WARN1 for specified commnicator.
template<typename... T_LOGOG_ARGS> void warn1(MPI_Comm comm, T_LOGOG_ARGS const& ... args)
{
    int rank;
    MPI_Comm_rank(comm, &rank);
    if(rank == 0)
        WARN1(args...);
};

/// Wrapped WARN2 for specified commnicator.
template<typename... T_LOGOG_ARGS> void warn2(MPI_Comm comm, T_LOGOG_ARGS const& ... args)
{
    int rank;
    MPI_Comm_rank(comm, &rank);
    if(rank == 0)
        WARN2(args...);
};

/// Wrapped WARN3 for specified commnicator.
template<typename... T_LOGOG_ARGS> void warn3(MPI_Comm comm, T_LOGOG_ARGS const& ... args)
{
    int rank;
    MPI_Comm_rank(comm, &rank);
    if(rank == 0)
        WARN3(args...);
};

/// Wrapped DBUG for specified commnicator.
template<typename... T_LOGOG_ARGS> void dbug(MPI_Comm comm, T_LOGOG_ARGS const& ... args)
{
    int rank;
    MPI_Comm_rank(comm, &rank);
    if(rank == 0)
        DBUG(args...);
};

/// Wrapped ERR for specified commnicator.
template<typename... T_LOGOG_ARGS> void err(MPI_Comm comm, T_LOGOG_ARGS const& ... args)
{
    int rank;
    MPI_Comm_rank(comm, &rank);
    if(rank == 0)
        ERR(args...);
};

/// Wrapped ALERT for specified commnicator.
template<typename... T_LOGOG_ARGS> void alert(MPI_Comm comm, T_LOGOG_ARGS const& ... args)
{
    int rank;
    MPI_Comm_rank(comm, &rank);
    if(rank == 0)
        ALERT(args...);
};

/// Wrapped CRITICAL for specified commnicator.
template<typename... T_LOGOG_ARGS> void critical(MPI_Comm comm, T_LOGOG_ARGS const& ... args)
{
    int rank;
    MPI_Comm_rank(comm, &rank);
    if(rank == 0)
        CRITICAL(args...);
}

/// Wrapped EMERGENCY for specified commnicator.
template<typename... T_LOGOG_ARGS> void emergency(MPI_Comm comm, T_LOGOG_ARGS const& ... args)
{
    int rank;
    MPI_Comm_rank(comm, &rank);
    if(rank == 0)
        EMERGENCY(args...);
};
#endif //USE_MPI

} // end namespace BaseLib

#endif

