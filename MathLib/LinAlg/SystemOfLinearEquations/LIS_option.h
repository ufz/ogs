/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 * \file LIS_option.h
 *
 * Created on 2012-06-25 by Norihiro Watanabe
 */

#ifndef LIS_OPTION_H_
#define LIS_OPTION_H_

#include <string>

#include "BaseLib/CodingTools.h"

/*
LIS matrix solver options
CG -i {cg|1}
BiCG -i {bicg|2}
CGS -i {cgs|3}
BiCGSTAB -i {bicgstab|4}
BiCGSTAB(l) -i {bicgstabl|5} -ell [2] Value for l
GPBiCG -i {gpbicg|6}
TFQMR -i {tfqmr|7}
Orthomin(m) -i {orthomin|8} -restart [40] Value for Restart m
GMRES(m) -i {gmres|9} -restart [40] Value for Restart m
Jacobi -i {jacobi|10}
Gauss-Seidel -i {gs|11}
SOR -i {sor|12} -omega [1.9] Value for Relaxation Coefficient  (0 <  < 2)
BiCGSafe -i {bicgsafe|13}
CR -i {cr|14}
BiCR -i {bicr|15}
CRS -i {crs|16}
BiCRSTAB -i {bicrstab|17}
GPBiCR -i {gpbicr|18}
BiCRSafe -i {bicrsafe|19}
FGMRES(m) -i {fgmres|20} -restart [40] Value for Restart m
IDR(s) -i {idrs|21} -restart [40] Value for Restart s

Preconditioner Option Auxiliary Option
None -p {none|0}
Jacobi -p {jacobi|1}
ILU(k) -p {ilu|2} -ilu_fill [0] Fill level k
SSOR -p {ssor|3} -ssor_w [1.0] Relaxation Coefficient  (0 <  < 2)
Hybrid -p {hybrid|4} -hybrid_i [sor] Iterative method
-hybrid_maxiter [25] Maximum number of iterations
-hybrid_tol [1.0e-3] Convergence criteria
-hybrid_w [1.5] Relaxation Coefficient  for
the SOR method (0 <  < 2)
-hybrid_ell [2] Value for l of the BiCGSTAB(l) method
-hybrid_restart [40] Restart values for GMRES and Orthomin
I+S -p {is|5} -is_alpha [1.0] Parameter ?for preconditioner
of a I + ?(m) type
-is_m [3] Parameter m for preconditioner
of a I + ?(m) type
SAINV -p {sainv|6} -sainv_drop [0.05] Drop criteria
SA-AMG -p {saamg|7} -saamg_unsym [false] Selection of asymmetric version
Crout ILU -p {iluc|8} -iluc_drop [0.05] Drop criteria
-iluc_rate [5.0] Ratio of Maximum fill-in
ILUT -p {ilut|9} -ilut_drop [0.05] Drop criteria
-ilut_rate [5.0] Ratio of Maximum fill-in
additive Schwarz -adds true -adds_iter [1] Number of iterations
*/

namespace MathLib
{

/**
 * \brief Option for Lis solver
 */
struct LIS_option
{
    enum SolverType
    {
        CG = 1,
        BiCG = 2,
        CGS = 3,
        BiCGSTAB = 4,
        BiCGSTABl = 5,
        GPBiCG = 6,
        TFQMR = 7,
        Orthomin = 8,
        GMRES = 9,
        INVALID = 0
    };

    enum PreconType
    {
        NONE = 0,
        JACOBI = 1,
        ILU = 2
    };

    int ls_method;
    int ls_precond;
    std::string ls_extra_arg;
    long ls_max_iterations;
    double ls_error_tolerance;

    LIS_option()
    {
        ls_method = CG;
        ls_precond = NONE;
        ls_max_iterations = 500;
        ls_error_tolerance = 1.e-6;
    }

    SolverType getSolverType(const std::string &str)
    {
        if (str.compare("CG")==0)
            return CG;
        if (str.compare("BICG")==0)
            return BiCG;
        if (str.compare("CGS")==0)
            return CGS;
        if (str.compare("BICGSTAB")==0)
            return BiCGSTAB;
        if (str.compare("BiCGSTABl")==0)
            return BiCGSTABl;
        if (str.compare("GPBiCG")==0)
            return GPBiCG;
        if (str.compare("TFQMR")==0)
            return TFQMR;
        if (str.compare("Orthomin")==0)
            return Orthomin;

        return INVALID;
    }

    PreconType getPreconType(const std::string &str)
    {
        if (str.compare("NONE")==0)
            return NONE;
        if (str.compare("JACOBI")==0)
            return JACOBI;
        if (str.compare("ILU")==0)
            return ILU;

        return NONE;
    }
};

}
#endif //LIS_OPTION_H_
