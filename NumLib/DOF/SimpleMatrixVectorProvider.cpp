/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "SimpleMatrixVectorProvider.h"

#include <cassert>
#include <logog/include/logog.hpp>

#include "BaseLib/Error.h"
#include "MathLib/LinAlg/LinAlg.h"
#include "MathLib/LinAlg/MatrixVectorTraits.h"

namespace LinAlg = MathLib::LinAlg;

namespace detail
{

template<typename MatVec>
MatVec*
transfer(std::map<std::size_t, MatVec*>& from_unused,
         std::map<MatVec*, std::size_t>& to_used,
         typename std::map<std::size_t, MatVec*>::iterator it)
{
    auto const id = it->first;
    auto& ptr = it->second;

    auto res = to_used.emplace(ptr, id);
    assert(res.second && "Emplacement failed.");
    from_unused.erase(it);
    return res.first->first;
}

template<typename MatVec>
void
transfer(std::map<MatVec*, std::size_t>& from_used,
         std::map<std::size_t, MatVec*>& to_unused,
         typename std::map<MatVec*, std::size_t>::iterator it)
{
    auto& ptr = it->first;
    auto const id = it->second;

    auto res = to_unused.emplace(id, ptr);
    assert(res.second && "Emplacement failed.");
    (void) res; // res unused if NDEBUG
    from_used.erase(it);
}

} // detail


namespace NumLib
{

template<bool do_search, typename MatVec, typename... Args>
std::pair<MatVec*, bool>
SimpleMatrixVectorProvider::
get_(std::size_t& id,
     std::map<std::size_t, MatVec*>& unused_map,
     std::map<MatVec*, std::size_t>& used_map,
     Args&&... args)
{
    if (id >= _next_id) {
        OGS_FATAL("An obviously uninitialized id argument has been passed."
            " This might not be a serious error for the current implementation,"
            " but it might become one in the future."
            " Hence, I will abort now.");
    }

    if (do_search)
    {
        auto it = unused_map.find(id);
        if (it != unused_map.end()) // unused matrix/vector found
            return { ::detail::transfer(unused_map, used_map, it), false };
    }

    // not searched or not found, so create a new one
    id = _next_id++;
    auto res = used_map.emplace(
        MathLib::MatrixVectorTraits<MatVec>::newInstance(std::forward<Args>(args)...).release(),
        id);
    assert(res.second && "Emplacement failed.");
    return { res.first->first, true };
}

template<bool do_search, typename... Args>
std::pair<GlobalMatrix*, bool>
SimpleMatrixVectorProvider::
getMatrix_(std::size_t& id, Args&&... args)
{
    return get_<do_search>(id, _unused_matrices, _used_matrices, std::forward<Args>(args)...);
}


GlobalMatrix&
SimpleMatrixVectorProvider::
getMatrix()
{
    std::size_t id = 0u;
    return *getMatrix_<false>(id).first;
}

GlobalMatrix&
SimpleMatrixVectorProvider::
getMatrix(std::size_t& id)
{
    return *getMatrix_<true>(id).first;
}

GlobalMatrix&
SimpleMatrixVectorProvider::
getMatrix(MathLib::MatrixSpecifications const& ms)
{
    std::size_t id = 0u;
    return *getMatrix_<false>(id, ms).first;
    // TODO assert that the returned object always is of the right size
}

GlobalMatrix&
SimpleMatrixVectorProvider::
getMatrix(MathLib::MatrixSpecifications const& ms, std::size_t& id)
{
    return *getMatrix_<true>(id, ms).first;
    // TODO assert that the returned object always is of the right size
}

GlobalMatrix&
SimpleMatrixVectorProvider::
getMatrix(GlobalMatrix const& A)
{
    std::size_t id = 0u;
    auto const& res = getMatrix_<false>(id, A);
    if (!res.second) // no new object has been created
        LinAlg::copy(A, *res.first);
    return *res.first;
}

GlobalMatrix&
SimpleMatrixVectorProvider::
getMatrix(GlobalMatrix const& A, std::size_t& id)
{
    auto const& res = getMatrix_<true>(id, A);
    if (!res.second) // no new object has been created
        LinAlg::copy(A, *res.first);
    return *res.first;
}

void
SimpleMatrixVectorProvider::
releaseMatrix(GlobalMatrix const& A)
{
    auto it = _used_matrices.find(const_cast<GlobalMatrix*>(&A));
    if (it == _used_matrices.end()) {
        OGS_FATAL("The given matrix has not been found. Cannot release it. Aborting.");
    } else {
        ::detail::transfer(_used_matrices, _unused_matrices, it);
    }
}

template<bool do_search, typename... Args>
std::pair<GlobalVector*, bool>
SimpleMatrixVectorProvider::
getVector_(std::size_t& id, Args&&... args)
{
    return get_<do_search>(id, _unused_vectors, _used_vectors, std::forward<Args>(args)...);
}


GlobalVector&
SimpleMatrixVectorProvider::
getVector()
{
    std::size_t id = 0u;
    return *getVector_<false>(id).first;
}

GlobalVector&
SimpleMatrixVectorProvider::
getVector(std::size_t& id)
{
    return *getVector_<true>(id).first;
}

GlobalVector&
SimpleMatrixVectorProvider::
getVector(MathLib::MatrixSpecifications const& ms)
{
    std::size_t id = 0u;
    return *getVector_<false>(id, ms).first;
    // TODO assert that the returned object always is of the right size
}

GlobalVector&
SimpleMatrixVectorProvider::
getVector(MathLib::MatrixSpecifications const& ms, std::size_t& id)
{
    return *getVector_<true>(id, ms).first;
    // TODO assert that the returned object always is of the right size
}

GlobalVector&
SimpleMatrixVectorProvider::
getVector(GlobalVector const& x)
{
    std::size_t id = 0u;
    auto const& res = getVector_<false>(id, x);
    if (!res.second) // no new object has been created
        LinAlg::copy(x, *res.first);
    return *res.first;
}

GlobalVector&
SimpleMatrixVectorProvider::
getVector(GlobalVector const& x, std::size_t& id)
{
    auto const& res = getVector_<true>(id, x);
    if (!res.second) // no new object has been created
        LinAlg::copy(x, *res.first);
    return *res.first;
}

void
SimpleMatrixVectorProvider::
releaseVector(GlobalVector const& x)
{
    auto it = _used_vectors.find(const_cast<GlobalVector*>(&x));
    if (it == _used_vectors.end()) {
        OGS_FATAL("The given vector has not been found. Cannot release it. Aborting.");
    } else {
        ::detail::transfer(_used_vectors, _unused_vectors, it);
    }
}

SimpleMatrixVectorProvider::
~SimpleMatrixVectorProvider()
{
    if ((!_used_matrices.empty()) || (!_used_vectors.empty())) {
        WARN("There are still some matrices and vectors in use."
             " This might be an indicator of a possible waste of memory.");
    }

    for (auto& id_ptr : _unused_matrices)
        delete id_ptr.second;

    for (auto& ptr_id : _used_matrices)
        delete ptr_id.first;

    for (auto& id_ptr : _unused_vectors)
        delete id_ptr.second;

    for (auto& ptr_id : _used_vectors)
        delete ptr_id.first;
}

} // MathLib
