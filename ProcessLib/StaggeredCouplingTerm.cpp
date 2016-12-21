/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
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
    return StaggeredCouplingTerm(coupled_pcsss, coupled_xs);
}

const LocalCouplingTerm createVoidLocalCouplingTerm()
{
    std::map<ProcessType, Process const&> coupled_pcsss;
    std::map<ProcessType, const std::vector<double>> local_coupled_xs;
    return LocalCouplingTerm(coupled_pcsss, std::move(local_coupled_xs));
}

} // end of ProcessLib

