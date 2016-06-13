/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <memory>
#include <gtest/gtest.h>

#include "NumLib/DOF/MeshComponentMap.h"

#include "MeshLib/MeshGenerators/MeshGenerator.h"
#include "MeshLib/MeshSubsets.h"

#include "NumLib/Assembler/VectorMatrixBuilder.h"

template <typename Builder>
class NumLibVectorMatrixBuilder : public ::testing::Test
{
    public:
    typedef MeshLib::MeshItemType MeshItemType;
    typedef MeshLib::Location Location;
    typedef NumLib::MeshComponentMap MeshComponentMap;

    typedef typename Builder::VectorType VectorType;
    typedef typename Builder::MatrixType MatrixType;

    public:
    NumLibVectorMatrixBuilder()
        : mesh(nullptr), nodesSubset(nullptr), cmap(nullptr)
    {
        mesh = MeshLib::MeshGenerator::generateLineMesh(1.0, mesh_size);
        nodesSubset = new MeshLib::MeshSubset(*mesh, &mesh->getNodes());

        // Add two components both based on the same nodesSubset.
        components.emplace_back(new MeshLib::MeshSubsets{nodesSubset});
        components.emplace_back(new MeshLib::MeshSubsets{nodesSubset});

        cmap = new MeshComponentMap(components,
            NumLib::ComponentOrder::BY_COMPONENT);
    }

    ~NumLibVectorMatrixBuilder()
    {
        delete cmap;
        delete nodesSubset;
        delete mesh;
    }

    static std::size_t const mesh_size = 9;
    MeshLib::Mesh const* mesh;
    MeshLib::MeshSubset const* nodesSubset;

    std::vector<std::unique_ptr<MeshLib::MeshSubsets>> components;
    MeshComponentMap const* cmap;
};

TYPED_TEST_CASE_P(NumLibVectorMatrixBuilder);

#ifndef USE_PETSC
TYPED_TEST_P(NumLibVectorMatrixBuilder, createVector)
#else
TYPED_TEST_P(NumLibVectorMatrixBuilder, DISABLED_createVector)
#endif
{
    typedef typename TestFixture::VectorType V;
    typedef TypeParam Builder;
    V* v = Builder::createVector(this->cmap->dofSizeWithGhosts());

    ASSERT_TRUE(v != nullptr);
    ASSERT_EQ(this->cmap->dofSizeWithGhosts(), v->size());

    delete v;
}

#ifndef USE_PETSC
TYPED_TEST_P(NumLibVectorMatrixBuilder, createMatrix)
#else
TYPED_TEST_P(NumLibVectorMatrixBuilder, DISABLED_createMatrix)
#endif
{
    typedef typename TestFixture::MatrixType M;
    typedef TypeParam Builder;
    M* m = Builder::createMatrix(this->cmap->dofSizeWithGhosts());

    ASSERT_TRUE(m != nullptr);
    ASSERT_EQ(this->cmap->dofSizeWithGhosts(), m->getNumberOfRows());
    ASSERT_EQ(this->cmap->dofSizeWithGhosts(), m->getNumberOfColumns());

    delete m;
}

#ifndef USE_PETSC
REGISTER_TYPED_TEST_CASE_P(NumLibVectorMatrixBuilder,
    createVector, createMatrix);
#else
REGISTER_TYPED_TEST_CASE_P(NumLibVectorMatrixBuilder,
    DISABLED_createVector, DISABLED_createMatrix);
#endif

#ifdef USE_LIS
#include "MathLib/LinAlg/Lis/LisVector.h"
#include "MathLib/LinAlg/Lis/LisMatrix.h"
#endif  // USE_LIS

#ifdef USE_PETSC
#include "MathLib/LinAlg/PETSc/PETScVector.h"
#include "MathLib/LinAlg/PETSc/PETScMatrix.h"
#endif  // USE_PETSC

#ifdef OGS_USE_EIGEN
#include "MathLib/LinAlg/Eigen/EigenVector.h"
#include "MathLib/LinAlg/Eigen/EigenMatrix.h"
#endif  // OGS_USE_EIGEN

typedef ::testing::Types
    <
      NumLib::VectorMatrixBuilder<
      MathLib::EigenMatrix, MathLib::EigenVector>
#ifdef USE_LIS
    , NumLib::VectorMatrixBuilder<
        MathLib::LisMatrix, MathLib::LisVector>
#endif  // USE_LIS
#ifdef USE_PETSC
    , NumLib::VectorMatrixBuilder<
        MathLib::PETScMatrix, MathLib::PETScVector>
#endif  // USE_PETSC
    > TestTypes;

INSTANTIATE_TYPED_TEST_CASE_P(templated, NumLibVectorMatrixBuilder,
    TestTypes);
