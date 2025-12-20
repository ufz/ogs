// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "PETScTools.h"

namespace MathLib
{
void applyKnownSolution(PETScMatrix& A, PETScVector& b, PETScVector& x,
                        const std::vector<PetscInt>& vec_knownX_id,
                        const std::vector<PetscScalar>& vec_knownX_x,
                        DirichletBCApplicationMode const /*mode*/)
{
    A.finalizeAssembly();

#ifndef NDEBUG
    if (std::any_of(vec_knownX_id.begin(), vec_knownX_id.end(),
                    [](PetscInt const i) { return i < 0; }))
    {
        OGS_FATAL(
            "Found negative indices in the vector of Dirichlet boundary "
            "conditions.");
    }
#endif  // NDEBUG

    A.setRowsColumnsZero(vec_knownX_id);
    A.finalizeAssembly();

    x.finalizeAssembly();
    b.finalizeAssembly();
    if (vec_knownX_id.size() > 0)
    {
        x.set(vec_knownX_id, vec_knownX_x);
        b.set(vec_knownX_id, vec_knownX_x);
    }

    x.finalizeAssembly();
    b.finalizeAssembly();
}

}  // end of namespace MathLib
