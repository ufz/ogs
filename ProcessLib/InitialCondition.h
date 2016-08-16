/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESS_LIB_INITIAL_CONDITION_H_
#define PROCESS_LIB_INITIAL_CONDITION_H_

#include "ProcessLib/Parameter/Parameter.h"

namespace ProcessLib
{
// TODO document
class InitialCondition
{
public:
    InitialCondition(Parameter<double> const& param) : _param(param) {}

    std::vector<double> const& getTuple(
            double const t, SpatialPosition const& pos) const
    {
        return _param.getTuple(t, pos);
    }

private:
    Parameter<double> const& _param;
};

}  // namespace ProcessLib

#endif  // PROCESS_LIB_INITIAL_CONDITION_H_
