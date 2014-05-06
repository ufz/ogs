/*!
   \file  PETScPC_KSP_GMRES_Option.cpp
   \brief Define the configuration data for the PETSc GMRES linear solver.

   \author Wenqing Wang
   \date 04-2014

   \copyright
    Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
               Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license
*/

#include "PETScPC_KSP_GMRES_Option.h"

namespace MathLib
{

using boost::property_tree::ptree;

PETScPC_KSP_GMRES_Option::
PETScPC_KSP_GMRES_Option(const boost::property_tree::ptree &option)
    : restart_number_gmres(30), is_modified_gram_schmidt_gmres(false),
      refine_type_gmres(KSP_GMRES_CGS_REFINE_NEVER)
{
    auto val = option.get_optional<double>("restart_number");
    restart_number_gmres = *val;

    boost::optional<bool> bool_vals = option.get_optional<bool>("is_modified_ram_schmidt_orthog");
    is_modified_gram_schmidt_gmres = *bool_vals;

    auto refine_type = option.get_optional<int>("refine_type");
    switch(*refine_type)
    {
        case 0:
            refine_type_gmres = KSP_GMRES_CGS_REFINE_NEVER;
            break;
        case 1:
            refine_type_gmres = KSP_GMRES_CGS_REFINE_IFNEEDED;
            break;
        case 2:
            refine_type_gmres = KSP_GMRES_CGS_REFINE_ALWAYS;
            break;
        default:
            refine_type_gmres = KSP_GMRES_CGS_REFINE_NEVER;
            break;
    }
}

/// Set Chebyshev option
void PETScPC_KSP_GMRES_Option::setOption(KSP &solver)
{
    KSPGMRESSetRestart(solver, restart_number_gmres);

    if(is_modified_gram_schmidt_gmres)
    {
        KSPGMRESSetOrthogonalization(solver, KSPGMRESClassicalGramSchmidtOrthogonalization);
    }

    KSPGMRESSetCGSRefinementType(solver, refine_type_gmres);
}

} // end namespace

