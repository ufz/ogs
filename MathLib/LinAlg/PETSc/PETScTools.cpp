/*!
   \file  PETScTools.cpp
   \brief Definition of a function related to PETSc solver interface to assign
         the Dirichlet boundary conditions.

   \author Wenqing Wang
   \version
   \date Nov 2011 - Sep 2013

   \copyright
    Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
               Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license
*/

#include "PETScTools.h"

#include "PETScMatrix.h"
#include "PETScVector.h"

#include "MeshLib/NodePartitionedMesh.h"
#include "MathLib/LinAlg/PETSc/PETScMatrixOption.h"

namespace MathLib
{

void applyKnownSolution(PETScMatrix &A, PETScVector &b, PETScVector &x,
                        const std::vector<PetscInt>  &vec_knownX_id,
                        const std::vector<PetscScalar> &vec_knownX_x)
{
    A.finalizeAssembly();

    A.setRowsColumnsZero(vec_knownX_id);
    A.finalizeAssembly();

    x.finalizeAssembly();
    b.finalizeAssembly();
    if(vec_knownX_id.size() > 0)
    {
        x.set(vec_knownX_id, vec_knownX_x);
        b.set(vec_knownX_id, vec_knownX_x);
    }

    x.finalizeAssembly();
    b.finalizeAssembly();
}

PETScMatrix* PETscMatrixAndNodeAdjacencyTableBuilder
::createMatrixAndNodeAdjacencyTable(const PetscInt dim, const MeshLib::Mesh& mesh,
                               MeshLib::NodeAdjacencyTable& /*node_adjacency_table*/)
{
    MathLib::PETScMatrixOption mat_opt;
    const MeshLib::NodePartitionedMesh &pmesh
                = static_cast<const MeshLib::NodePartitionedMesh&>(mesh);
    mat_opt.d_nz = pmesh.getMaximumNConnectedNodesToNode();
    mat_opt.o_nz = mat_opt.d_nz;

    return new PETScMatrix(dim, mat_opt);
}

} // end of namespace MathLib


