/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

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
class SimpleMatrixVectorProvider final
        : public MatrixProvider
        , public VectorProvider
{
public:
    SimpleMatrixVectorProvider() = default;

    // no copies
    SimpleMatrixVectorProvider(SimpleMatrixVectorProvider const&) = delete;
    SimpleMatrixVectorProvider& operator=(SimpleMatrixVectorProvider const&) = delete;

    GlobalVector& getVector() override;
    GlobalVector& getVector(std::size_t& id) override;

    GlobalVector& getVector(GlobalVector const& x) override;
    GlobalVector& getVector(GlobalVector const& x, std::size_t& id) override;

    GlobalVector& getVector(MathLib::MatrixSpecifications const& ms) override;
    GlobalVector& getVector(MathLib::MatrixSpecifications const& ms, std::size_t& id) override;

    void releaseVector(GlobalVector const& x) override;

    GlobalMatrix& getMatrix() override;
    GlobalMatrix& getMatrix(std::size_t& id) override;

    GlobalMatrix& getMatrix(GlobalMatrix const& A) override;
    GlobalMatrix& getMatrix(GlobalMatrix const& A, std::size_t& id) override;

    GlobalMatrix& getMatrix(MathLib::MatrixSpecifications const& ms) override;
    GlobalMatrix& getMatrix(MathLib::MatrixSpecifications const& ms, std::size_t& id) override;

    void releaseMatrix(GlobalMatrix const& A) override;

    ~SimpleMatrixVectorProvider() override;

private:
    template<bool do_search, typename... Args>
    std::pair<GlobalMatrix*, bool> getMatrix_(std::size_t& id, Args&&... args);

    template<bool do_search, typename... Args>
    std::pair<GlobalVector*, bool> getVector_(std::size_t& id, Args&&... args);

    // returns a pair with the pointer to the matrix/vector and
    // a boolean indicating if a new object has been built (then true else false)
    template<bool do_search, typename MatVec, typename... Args>
    std::pair<MatVec*, bool>
    get_(std::size_t& id,
         std::map<std::size_t, MatVec*>& unused_map,
         std::map<MatVec*, std::size_t>& used_map,
         Args&&... args);

    std::size_t _next_id = 1;

    std::map<std::size_t, GlobalMatrix*> _unused_matrices;
    std::map<GlobalMatrix*, std::size_t> _used_matrices;

    std::map<std::size_t, GlobalVector*> _unused_vectors;
    std::map<GlobalVector*, std::size_t> _used_vectors;
};



} // namespace NumLib
