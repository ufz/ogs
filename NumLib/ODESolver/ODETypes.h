// TODO: move that file somewhere else

#pragma once

#include <initializer_list>
#include <cassert>


// Eigen dense matrix ///////////////////////
// always enabled

inline void setMatrix(Eigen::MatrixXd& m,
                      Eigen::MatrixXd::Index const rows, Eigen::MatrixXd::Index const cols,
                      std::initializer_list<double> values)
{
    using IndexType = Eigen::MatrixXd::Index;
    assert((IndexType) values.size() == rows*cols);

    auto it = values.begin();
    for (IndexType r=0; r<rows; ++r) {
        for (IndexType c=0; c<cols; ++c) {
            m(r, c) = *(it++);
        }
    }
}

inline void setMatrix(Eigen::MatrixXd& m, Eigen::MatrixXd const& tmp)
{
    m = tmp;
}

inline void addToMatrix(Eigen::MatrixXd& m,
                        Eigen::MatrixXd::Index const rows, Eigen::MatrixXd::Index const cols,
                        std::initializer_list<double> values)
{
    using IndexType = Eigen::MatrixXd::Index;
    assert((IndexType) values.size() == rows*cols);

    auto it = values.begin();
    for (IndexType r=0; r<rows; ++r) {
        for (IndexType c=0; c<cols; ++c) {
            m(r, c) += *(it++);
        }
    }
}

namespace NumLib
{
inline void oneShotLinearSolve(Eigen::MatrixXd& A, Eigen::VectorXd& rhs, Eigen::VectorXd& x)
{
    // Eigen::SparseLU<Eigen::MatrixXd> slv;
    Eigen::FullPivLU<Eigen::MatrixXd> slv;
    slv.compute(A);
    x = slv.solve(rhs);
}

inline double norm(Eigen::VectorXd const& x) { return x.norm(); }

}

inline void setVector(Eigen::VectorXd& v, std::initializer_list<double> values)
{
    assert((std::size_t) v.size() == values.size());
    auto it = values.begin();
    for (std::size_t i=0; i<values.size(); ++i) v[i] = *(it++);
}








/*
#i fdef USE_EIGEN_PLAIN
#include <Eigen/SparseCore>
// #include <Eigen/SparseLU>
#include <Eigen/LU>

using IndexType = int;

namespace NumLib
{
inline void oneShotLinearSolve(Eigen::MatrixXd& A, Eigen::VectorXd& rhs, Eigen::VectorXd& x)
{
    // Eigen::SparseLU<Eigen::MatrixXd> slv;
    Eigen::FullPivLU<Eigen::MatrixXd> slv;
    slv.compute(A);
    x = slv.solve(rhs);
}
}

*/


#ifdef OGS_USE_EIGEN

#include "MathLib/LinAlg/Eigen/EigenVector.h"
#include "MathLib/LinAlg/Eigen/EigenMatrix.h"
#include "MathLib/LinAlg/Eigen/EigenLinearSolver.h"

using IndexType = int;

namespace NumLib
{
inline void oneShotLinearSolve(MathLib::EigenMatrix& A, MathLib::EigenVector& rhs, MathLib::EigenVector& x)
{
    MathLib::EigenLinearSolver slv(A);
    slv.solve(rhs, x);
}
}


inline void setVector(MathLib::EigenVector& v,
                      std::initializer_list<double> values)
{
    setVector(v.getRawVector(), values);
}


inline void setMatrix(MathLib::EigenMatrix& m,
                      IndexType const rows, IndexType const cols,
                      std::initializer_list<double> values)
{
    assert((IndexType) values.size() == rows*cols);
    Eigen::MatrixXd tmp(rows, cols);

    auto it = values.begin();
    for (IndexType r=0; r<rows; ++r) {
        for (IndexType c=0; c<cols; ++c) {
            tmp(r, c) = *(it++);
        }
    }

    m.getRawMatrix() = tmp.sparseView();
}

inline void setMatrix(MathLib::EigenMatrix& m, Eigen::MatrixXd const& tmp)
{
    m.getRawMatrix() = tmp.sparseView();
}

inline void addToMatrix(MathLib::EigenMatrix& m,
                        IndexType const rows, IndexType const cols,
                        std::initializer_list<double> values)
{
    assert((IndexType) values.size() == rows*cols);
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> tmp(rows, cols);

    auto it = values.begin();
    for (IndexType r=0; r<rows; ++r) {
        for (IndexType c=0; c<cols; ++c) {
            tmp(r, c) = *(it++);
        }
    }

    m.getRawMatrix() += tmp.sparseView();
}



#endif // OGS_USE_EIGEN
