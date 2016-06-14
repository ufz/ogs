/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "LinearSolver.h"

namespace
{
//! Describes the default linear solver type to be used together with the given
//! \c Matrix and \c Vector types.
template<typename Matrix, typename Vector>
struct DefaultLinearSolver;
}

namespace MathLib
{
template <typename Matrix, typename Vector, typename Solver>
std::unique_ptr<LinearSolver<Matrix, Vector>> createLinearSolver(
    BaseLib::ConfigTree const* const config)
{
    using AbstractLS = LinearSolver<Matrix, Vector>;
    return std::unique_ptr<AbstractLS>{new Solver{"", config}};
}

template <typename Matrix, typename Vector>
std::unique_ptr<LinearSolver<Matrix, Vector>> createLinearSolver(
    BaseLib::ConfigTree const* const config)
{
    using Solver = typename DefaultLinearSolver<Matrix, Vector>::Solver;
    return createLinearSolver<Matrix, Vector, Solver>(config);
}

} // MathLib


// specializes the DefaultLinearSolver<> class
#define DEFAULT_LINEAR_SOLVER(MAT, VEC, SLV) \
    namespace { \
    template<> struct DefaultLinearSolver<MAT, VEC> { \
        using Solver = SLV; \
    };}

// does explicit template instantiation of createLinearSolver<>()
#define CREATE_LINEAR_SOLVER(MAT, VEC, SLV) \
    namespace MathLib { \
    template std::unique_ptr<LinearSolver<MAT, VEC> > \
    createLinearSolver<MAT, VEC, SLV>(BaseLib::ConfigTree const*const); \
    template std::unique_ptr<LinearSolver<MAT, VEC> > \
    createLinearSolver<MAT, VEC>(BaseLib::ConfigTree const*const); \
    }

// both of the above
#define SPECIALIZE_CREATE_LINEAR_SOLVER(MAT, VEC, SLV) \
    DEFAULT_LINEAR_SOLVER(MAT, VEC, SLV) \
    CREATE_LINEAR_SOLVER(MAT, VEC, SLV)


#ifdef OGS_USE_EIGEN

#include <Eigen/QR>

namespace MathLib
{

//! Matches Eigen solvers to our interface.
class EigenDenseSolverWrapper
        : public LinearSolver<Eigen::MatrixXd, Eigen::VectorXd>
{
public:
    EigenDenseSolverWrapper(std::string const&, BaseLib::ConfigTree const*const) {}

    bool solve(Eigen::MatrixXd &A, Eigen::VectorXd &b, Eigen::VectorXd &x) override
    {
        _solver.compute(A);
        x = _solver.solve(b);
        return true; // TODO add checks
    }

private:
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> _solver;
};

} // namespace MathLib

SPECIALIZE_CREATE_LINEAR_SOLVER(Eigen::MatrixXd, Eigen::VectorXd,
                                MathLib::EigenDenseSolverWrapper)

#endif


#ifdef OGS_USE_LIS

#include "MathLib/LinAlg/EigenLis/EigenLisLinearSolver.h"

SPECIALIZE_CREATE_LINEAR_SOLVER(MathLib::EigenMatrix, MathLib::EigenVector,
                                MathLib::EigenLisLinearSolver)

#elif defined(USE_PETSC)

#include "MathLib/LinAlg/PETSc/PETScLinearSolver.h"

SPECIALIZE_CREATE_LINEAR_SOLVER(MathLib::PETScMatrix, MathLib::PETScVector,
                                MathLib::PETScLinearSolver)

#elif defined(OGS_USE_EIGEN)

#include "MathLib/LinAlg/Eigen/EigenLinearSolver.h"

SPECIALIZE_CREATE_LINEAR_SOLVER(MathLib::EigenMatrix, MathLib::EigenVector,
                                MathLib::EigenLinearSolver)

#endif
