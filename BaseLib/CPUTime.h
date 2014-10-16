/**
 * \file
 * \author Thomas Fischer
 * \author Wenqing Wang
 * \date   2012-05-10, 2014.10.10
 * \brief  Definition of the CPUTime class.
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
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

/// Record CPU time
class CPUTime
{
    public:
        /// Start the timer.
        void start()
        {
            _timer = - clock()/static_cast<double>(CLOCKS_PER_SEC);
        }

        /// Get the epalsed time after started.
        double elapsed()
        {
            return _timer + clock()/static_cast<double>(CLOCKS_PER_SEC);
        }
    private:
        double _timer = 0.;
};

} // end namespace BaseLib

#endif

