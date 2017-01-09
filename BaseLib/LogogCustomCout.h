/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef LOGOGCUSTOMCOUT_H_
#define LOGOGCUSTOMCOUT_H_

#include <ostream>

#include <logog/include/logog.hpp>

#ifdef USE_MPI
#include <mpi.h>
#endif

namespace BaseLib
{

/// Custom target for logog output
class LogogCustomCout : public logog::Target
{
public:
#ifdef USE_MPI
    /**
     * Constructor when MPI is involved
     *
     * @param all_rank_output_level  Minimum level to output messages from all MPI processes
     * @param mpi_comm               MPI communicator
     */
    LogogCustomCout(LOGOG_LEVEL_TYPE all_rank_output_level = LOGOG_LEVEL_INFO, MPI_Comm mpi_comm = MPI_COMM_WORLD)
    : _all_rank_output_level(all_rank_output_level), _is_rank0 (getRank(mpi_comm)==0)
    {}
#endif

    virtual int Receive( const logog::Topic &topic )
    {
#ifdef USE_MPI
        if (topic.Level() > _all_rank_output_level && !_is_rank0)
            return 0;
#endif
        return logog::Target::Receive(topic);
    }

    virtual int Output( const LOGOG_STRING &data )
    {
        LOGOG_COUT << (const LOGOG_CHAR *)data << std::flush;
        return 0;
    }

private:
#ifdef USE_MPI
    int getRank(MPI_Comm mpi_comm) const
    {
        int rank = 0;
        MPI_Comm_rank(mpi_comm, &rank);
        return rank;
    }

    const LOGOG_LEVEL_TYPE _all_rank_output_level;
    const bool _is_rank0;
#endif
};

#endif // LOGOGCUSTOMCOUT_H_

} // namespace BaseLib
