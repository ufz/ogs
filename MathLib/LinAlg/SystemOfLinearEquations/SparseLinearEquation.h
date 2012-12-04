/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 * \file SparseLinearEquation.h
 *
 * Created on 2012-06-25 by Norihiro Watanabe
 */

#ifndef SPARSELINEAREQUATION_H_
#define SPARSELINEAREQUATION_H_

#include <vector>
#include <map>

#include "BaseLib/CodingTools.h"

#include "MathLib/LinAlg/Dense/Matrix.h"
#include "MathLib/LinAlg/Sparse/CRSMatrix.h"
#include "AbstractCRSLinearEquation.h"

namespace MathLib
{

/**
 * \brief Linear equation using sparse solvers developed by TF
 */
class SparseLinearEquation : public AbstractCRSLinearEquation<unsigned>
{
public:
    enum SolverType
    {
        INVALID,
        SolverCG,
        SolverBiCGStab,
        SolverGMRes
    };

    enum PreconditionerType
    {
        NONE,
        DIAG,
        DIAGSCALE
    };

    struct SpLinearOptions
    {
        SolverType solver_type;
        PreconditionerType precon_type;
        double error_tolerance;
        size_t max_iteration_step;

        SpLinearOptions()
        {
            solver_type = SolverCG;
            precon_type = NONE;
            error_tolerance = 1.e-10;
            max_iteration_step = 500;
        };


        SolverType getSolverType(const std::string &str)
        {
            if (str.compare("CG")==0)
                return SolverCG;
            if (str.compare("BICGSTAB")==0)
                return SolverBiCGStab;
            if (str.compare("GMRES")==0)
                return SolverGMRes;

            return INVALID;
        }
        PreconditionerType getPreconType(const std::string &str)
        {
            if (str.compare("NONE")==0)
                return NONE;
            if (str.compare("DIAG")==0)
                return DIAG;
            if (str.compare("DIAGSCALE")==0)
                return DIAGSCALE;

            return NONE;
        }
    };

    SparseLinearEquation(std::size_t length, RowMajorSparsity* sp) : AbstractCRSLinearEquation<unsigned>(length, sp) {};

    virtual ~SparseLinearEquation()
    {
    }

    void initialize() {};
    void finalize() {};

    void setOption(const BaseLib::Options &option);

    void setOption(const SpLinearOptions &option)
    {
        _option = option;
    }

    SpLinearOptions &getOption()
    {
        return _option;
    }


protected:
    void solveEqs(CRSMatrix<double, unsigned> *A, double *rhs, double *x);

private:
    SpLinearOptions _option;

    DISALLOW_COPY_AND_ASSIGN(SparseLinearEquation);
};

}

#endif //SPARSELINEAREQUATION_H_
