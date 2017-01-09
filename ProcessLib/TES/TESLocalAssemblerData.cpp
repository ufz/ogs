/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "TESLocalAssemblerData.h"
#include "TESReactionAdaptor.h"

namespace ProcessLib
{
namespace TES
{
TESLocalAssemblerData::TESLocalAssemblerData(AssemblyParams const& ap_,
                                             const unsigned num_int_pts,
                                             const unsigned dimension)
    : ap(ap_),
      solid_density(num_int_pts, ap_.initial_solid_density),
      reaction_rate(num_int_pts),
      velocity(dimension, std::vector<double>(num_int_pts)),
      reaction_adaptor(TESFEMReactionAdaptor::newInstance(*this)),
      solid_density_prev_ts(num_int_pts, ap_.initial_solid_density),
      reaction_rate_prev_ts(num_int_pts)
{
}

TESLocalAssemblerData::~TESLocalAssemblerData() = default;
}
}  // namespaces
