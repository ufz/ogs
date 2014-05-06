/*!
   \file  PETScPC_KSP_GMRES_Option.h
   \brief Define the configuration data for the PETSc GMRES linear solver.

   \author Wenqing Wang
   \date 04-2014

   \copyright
    Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
               Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license
*/
#ifndef PETSCPC_KSP_GMRES_OPTION_H_
#define PETSCPC_KSP_GMRES_OPTION_H_

#include <boost/property_tree/ptree.hpp>

#include <petscksp.h>

namespace MathLib
{
/*!
    PETSc KSP GMRES options

 */
struct PETScPC_KSP_GMRES_Option
{
    PETScPC_KSP_GMRES_Option(const boost::property_tree::ptree &option);

    /// Set GMRES option
    void setOption(KSP &solver);

    /// Overloaded assign operator
    void operator = (const PETScPC_KSP_GMRES_Option& opt)
    {
        restart_number_gmres = opt.restart_number_gmres;
        is_modified_gram_schmidt_gmres = opt.is_modified_gram_schmidt_gmres;
        refine_type_gmres = opt.refine_type_gmres;
    }

    ///  Restart number of GMRES.
    PetscInt restart_number_gmres;

    /// Flag for the modified Gram-Schmidt orthogonalization.
    bool is_modified_gram_schmidt_gmres;

    /*!
         \brief Refinement type for GMRES.
         This iterative refinement is used to improve the stability
         of orthogonalization.
         KSPGMRESCGSRefinementType is a enum type of PETSc defined as
         typedef enum {KSP_GMRES_CGS_REFINE_NEVER, KSP_GMRES_CGS_REFINE_IFNEEDED,
                       KSP_GMRES_CGS_REFINE_ALWAYS} KSPGMRESCGSRefinementType;
    */
    KSPGMRESCGSRefinementType refine_type_gmres;
};

} // end namespace
#endif

