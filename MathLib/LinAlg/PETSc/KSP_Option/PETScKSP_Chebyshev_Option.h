/*!
   \file  PETScKSP_Chebyshev_Option.h
   \brief Define the configuration data for the PETSc Chebyshev linear solver.

   \author Wenqing Wang
   \date 04-2014

   \copyright
    Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
               Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license
*/
#ifndef PETSCKSP_CHEBYSHEV_OPTION_H_
#define PETSCKSP_CHEBYSHEV_OPTION_H_

#include <boost/property_tree/ptree.hpp>

#include <petscksp.h>

namespace MathLib
{
/*!
    PETSc KSP Chebyshev options

 */
struct PETScKSP_Chebyshev_Option
{
    PETScKSP_Chebyshev_Option(const boost::property_tree::ptree &option);

    /// Set Chebyshev option
    void setOption(KSP &solver)
    {
        KSPChebyshevSetEigenvalues(solver, emax_chebyshev, emin_chebyshev);
    }

    /// Overloaded assign operator
    void operator = (const PETScKSP_Chebyshev_Option& opt)
    {
        emin_chebyshev = opt.emin_chebyshev;
        emax_chebyshev = opt.emax_chebyshev;
    }

    /// Smallest eignvalue for Chebyshev.
    PetscReal emin_chebyshev;
    /// maximum eignvalue for Chebyshev.
    PetscReal emax_chebyshev;
};

} // end namespace
#endif

