/**
 * \author Norihiro Watanabe
 * \date   2012-08-03
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include "BaseLib/Logging.h"

#include "NumLib/TimeStepping/TimeStep.h"


namespace
{

struct Dummy
{
    template <class T>
    void operator()(T &/*obj*/) {}
};

template <class T_TIME_STEPPING, class T = Dummy>
std::vector<double> timeStepping(T_TIME_STEPPING& algorithm,
                                 std::vector<int> const& number_iterations,
                                 T* obj = nullptr)
{
    std::vector<double> vec_t;
    vec_t.push_back(algorithm.begin());

    double const solution_error = 0;
    for (auto const& i : number_iterations)
    {
        if (!algorithm.next(solution_error, i))
        {
            break;
        }

        NumLib::TimeStep t = algorithm.getTimeStep();
        //INFO("t: n=%d,t=%g,dt=%g", t.steps(), t.current(), t.dt());
        if (obj)
        {
            (*obj)(algorithm);  // do something
        }
        if (algorithm.accepted()) {
            vec_t.push_back(t.current());
        } else {
            //INFO("*** rejected.");
        }
    }

    return vec_t;
}
} // namespace
