/**
 * \file
 * \author Norihiro Watanabe
 * \date   2012-06-25
 * \brief  Definition of the LisOption class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef LIS_OPTION_H_
#define LIS_OPTION_H_

#include <string>

namespace MathLib
{

/**
 * \brief Option for Lis solver
 */
struct LisOption
{
    /// Solver type
    enum class SolverType : int
    {
        INVALID = 0,
        CG = 1,
        BiCG = 2,
        CGS = 3,
        BiCGSTAB = 4,
        BiCGSTABl = 5,
        GPBiCG = 6,
        TFQMR = 7,
        Orthomin = 8,
        GMRES = 9,
        Jacobi = 10,
        GaussSeidel = 11,
        SOR = 12,
        BiCGSafe = 13,
        CR = 14,
        BiCR = 15,
        CRS = 16,
        BiCRSTAB = 17,
        GPBiCR = 18,
        BiCRSafe = 19,
        FGMRESm = 20,
        IDRs = 21,
        MINRES = 22
    };

    /// Preconditioner type
    enum class PreconType : int
    {
        NONE = 0,
        JACOBI = 1,
        ILU = 2,
        SSOR = 3,
        Hybrid = 4,
        IplusS = 5,
        SAINV = 6,
        SAAMG = 7,
        CroutILU = 8,
        ILUT = 9
    };

    /// Matrix type
    enum class MatrixType : int
    {
        CRS = 1,
        CCS = 2,
        MSR = 3,
        DIA = 4,
        ELL = 5,
        JDS = 6,
        BSR = 7,
        BSC = 8,
        VBR = 9,
        COO = 10,
        DNS = 11
    };

    /// Linear solver type
    SolverType solver_type;
    /// Preconditioner type
    PreconType precon_type;
    /// Matrix type
    MatrixType matrix_type;
    /// Maximum iteration count
    long max_iterations;
    /// Error tolerance
    double error_tolerance;
    /// Extra option
    std::string extra_arg;
    /// Arguments for solver and preconditioner. This variable is always preferred
    /// to other variables.
    std::string solver_precon_arg;

    /**
     * Constructor
     *
     * Default options are CG, no preconditioner, iteration count 500 and
     * tolerance 1e-10. Default matrix storage type is CRS.
     */
    LisOption();

    /// Destructor
    ~LisOption() {}

    /**
     * return a linear solver type from the solver name
     *
     * @param solver_name
     * @return a linear solver type
     *      If there is no solver type matched with the given name, INVALID
     *      is returned.
     */
    static SolverType getSolverType(const std::string &solver_name);

    /**
     * return a preconditioner type from the name
     *
     * @param precon_name
     * @return a preconditioner type
     *      If there is no preconditioner type matched with the given name, NONE
     *      is returned.
     */
    static PreconType getPreconType(const std::string &precon_name);

    /**
     * return a matrix type
     * @param matrix_name
     * @return
     *      If there is no matrix type matched with the given name, CRS
     *      is returned.
     */
    static MatrixType getMatrixType(const std::string &matrix_name);
};

}
#endif //LIS_OPTION_H_
