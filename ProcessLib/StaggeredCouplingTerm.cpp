/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file   StaggeredCouplingTerm.cpp
 *
 * Created on November 7, 2016, 12:14 PM
 */

#include "StaggeredCouplingTerm.h"
#include "Process.h"

namespace ProcessLib
{
const StaggeredCouplingTerm createVoidStaggeredCouplingTerm()
{
    std::map<ProcessType, Process const&> coupled_pcsss;
    std::map<ProcessType, GlobalVector const*> coupled_xs;
    const bool empty = true;
    return StaggeredCouplingTerm(coupled_pcsss, coupled_xs, empty);
}

} // end of ProcessLib

