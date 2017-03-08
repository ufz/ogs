/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "SmallDeformationProcess.h"

extern template class ProcessLib::LIE::SmallDeformation::SmallDeformationProcess<2>;
extern template class ProcessLib::LIE::SmallDeformation::
    SmallDeformationProcess<3>;
