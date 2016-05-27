/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef NUMLIB_TIMEDISCRETIZATION_BUILDER_H
#define NUMLIB_TIMEDISCRETIZATION_BUILDER_H

#include <memory>

#include "BaseLib/ConfigTree.h"
#include "TimeDiscretization.h"

namespace NumLib
{

template<typename Vector>
std::unique_ptr<TimeDiscretization<Vector> >
createTimeDiscretization(BaseLib::ConfigTree const& config)
{
    using T = std::unique_ptr<TimeDiscretization<Vector> >;

    //! \ogs_file_param{process__time_discretization__type}
    auto const type = config.getConfParam<std::string>("type");

    if (type == "BackwardEuler") {
        using ConcreteTD = BackwardEuler<Vector>;
        return T(new ConcreteTD);
    } else if (type == "ForwardEuler") {
        using ConcreteTD = ForwardEuler<Vector>;
        return T(new ConcreteTD);
    } else if (type == "CrankNicolson") {
        //! \ogs_file_param{process__time_discretization__CrankNicolson__theta}
        auto const theta = config.getConfParam<double>("theta");
        using ConcreteTD = CrankNicolson<Vector>;
        return T(new ConcreteTD(theta));
    } else if (type == "BackwardDifferentiationFormula") {
        //! \ogs_file_param{process__time_discretization__BackwardDifferentiationFormula__order}
        auto const order = config.getConfParam<unsigned>("order");
        using ConcreteTD = BackwardDifferentiationFormula<Vector>;
        return T(new ConcreteTD(order));
    } else {
        ERR("Unrecognized time discretization type `%s'", type.c_str());
        std::abort();
    }
}

}

#endif // NUMLIB_TIMEDISCRETIZATION_BUILDER_H
