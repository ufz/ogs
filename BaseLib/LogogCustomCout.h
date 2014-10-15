/**
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef LOGOGCUSTOMCOUT_H_
#define LOGOGCUSTOMCOUT_H_

#include "logog/include/logog.hpp"

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
	LogogCustomCout(MPI_Comm mpi_comm)
	{
		MPI_Comm_rank(mpi_comm, &_mpi_rank);
	}
#endif

	virtual int Output( const LOGOG_STRING &data )
	{
#ifdef USE_MPI
		if (_mpi_rank == 0)
#endif
		LOGOG_COUT << (const LOGOG_CHAR *)data;
		return 0;
	}

#ifdef USE_MPI
	int _mpi_rank;
#endif
};

#endif // LOGOGCUSTOMCOUT_H_

} // namespace BaseLib
