/**
 * \file
 * \author Thomas Fischer
 * \author Wenqing Wang
 * \date   2012-05-10, 2014.10.10
 * \brief  Definition of the CPUTime class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef CPUTIME_H
#define CPUTIME_H

#include <ctime>

namespace BaseLib
{

/// Count CPU time
class CPUTime
{
    public:
        /// Start the timer.
        void start()
        {
            _start_time = clock();
        }

        /// Get the elapsed time after started.
        double elapsed() const
        {
            return (clock() - _start_time)/static_cast<double>(CLOCKS_PER_SEC);
        }
    private:
        double _start_time = 0.;
};

} // end namespace BaseLib

#endif

