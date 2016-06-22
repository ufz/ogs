/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef NUMLIB_SIMPLE_MATRIX_PROVIDER_H
#define NUMLIB_SIMPLE_MATRIX_PROVIDER_H

#include <map>
#include <memory>

#include "MatrixProviderUser.h"

namespace NumLib
{

/*! Manages storage for matrices and vectors.
 *
 * This is a simple implementation of the MatrixProvider and VectorProvider interfaces.
 *
 * It is simple insofar it does not reuse released matrices/vectors, but keeps them in
 * memory until they are acquired again by the user.
 */
template<typename Matrix, typename Vector>
class SimpleMatrixVectorProvider final
        : public MatrixProvider<Matrix>
        , public VectorProvider<Vector>
{
public:
    SimpleMatrixVectorProvider() = default;

    // no copies
    SimpleMatrixVectorProvider(SimpleMatrixVectorProvider const&) = delete;
    SimpleMatrixVectorProvider& operator=(SimpleMatrixVectorProvider const&) = delete;

    Vector& getVector() override;
    Vector& getVector(std::size_t& id) override;

    Vector& getVector(Vector const& x) override;
    Vector& getVector(Vector const& x, std::size_t& id) override;

    Vector& getVector(MathLib::MatrixSpecifications const& ms) override;
    Vector& getVector(MathLib::MatrixSpecifications const& ms, std::size_t& id) override;

    void releaseVector(Vector const& x) override;

    Matrix& getMatrix() override;
    Matrix& getMatrix(std::size_t& id) override;

    Matrix& getMatrix(Matrix const& A) override;
    Matrix& getMatrix(Matrix const& A, std::size_t& id) override;

    Matrix& getMatrix(MathLib::MatrixSpecifications const& ms) override;
    Matrix& getMatrix(MathLib::MatrixSpecifications const& ms, std::size_t& id) override;

    void releaseMatrix(Matrix const& A) override;

    ~SimpleMatrixVectorProvider();

private:
    template<bool do_search, typename... Args>
    std::pair<Matrix*, bool> getMatrix_(std::size_t& id, Args&&... args);

    template<bool do_search, typename... Args>
    std::pair<Vector*, bool> getVector_(std::size_t& id, Args&&... args);

    // returns a pair with the pointer to the matrix/vector and
    // a boolean indicating if a new object has been built (then true else false)
    template<bool do_search, typename MatVec, typename... Args>
    std::pair<MatVec*, bool>
    get_(std::size_t& id,
         std::map<std::size_t, MatVec*>& unused_map,
         std::map<MatVec*, std::size_t>& used_map,
         Args&&... args);

    std::size_t _next_id = 1;

    std::map<std::size_t, Matrix*> _unused_matrices;
    std::map<Matrix*, std::size_t> _used_matrices;

    std::map<std::size_t, Vector*> _unused_vectors;
    std::map<Vector*, std::size_t> _used_vectors;
};



} // namespace NumLib

#include "SimpleMatrixVectorProvider-impl.h"

#endif // NUMLIB_SIMPLE_MATRIX_PROVIDER_H
