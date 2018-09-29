/**
 * \author Norihiro Watanabe
 * \date   2012-08-03
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <logog/include/logog.hpp>

#include "NumLib/TimeStepping/TimeStep.h"


namespace
{

struct Dummy
{
    template <class T>
    void operator()(T &/*obj*/) {}
};

template <class T_TIME_STEPPING, class T=Dummy>
std::vector<double> timeStepping(T_TIME_STEPPING &algorithm, T* obj=nullptr)
{
    std::vector<double> vec_t;
    vec_t.push_back(algorithm.begin());

    while (algorithm.next(0.0))
    {
        NumLib::TimeStep t = algorithm.getTimeStep();
        //INFO("t: n=%d,t=%g,dt=%g", t.steps(), t.current(), t.dt());
        if (obj)
            (*obj)(algorithm); // do something
        if (algorithm.accepted()) {
            vec_t.push_back(t.current());
        } else {
            //INFO("*** rejected.");
        }
    }

    return vec_t;
}
} // namespace
