// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <memory>

#include "MathLib/LinAlg/GlobalMatrixVectorTypes.h"
#include "MathLib/LinAlg/MatrixVectorTraits.h"

namespace ProcessLib
{
//! Stores assembled global M, K, b matrices/vectors for later reuse.
struct AssembledMatrixCache final
{
    //! Access the stored data.
    std::tuple<GlobalMatrix&, GlobalMatrix&, GlobalVector&> MKb() const
    {
        return std::tie(*M_, *K_, *b_);
    }

    //! Check if data has been stored.
    bool hasMKb() const { return M_ && K_ && b_; }

    //! Store data.
    void storeMKb(GlobalMatrix const& M,
                  GlobalMatrix const& K,
                  GlobalVector const& b)
    {
        M_ = MathLib::MatrixVectorTraits<GlobalMatrix>::newInstance(M);
        K_ = MathLib::MatrixVectorTraits<GlobalMatrix>::newInstance(K);
        b_ = MathLib::MatrixVectorTraits<GlobalVector>::newInstance(b);
    }

private:
    std::unique_ptr<GlobalMatrix> M_;
    std::unique_ptr<GlobalMatrix> K_;
    std::unique_ptr<GlobalVector> b_;
};
}  // namespace ProcessLib
