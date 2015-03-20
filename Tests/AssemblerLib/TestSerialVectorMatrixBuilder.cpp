/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <gtest/gtest.h>

#include "AssemblerLib/MeshComponentMap.h"

#include "MeshLib/MeshGenerators/MeshGenerator.h"
#include "MeshLib/MeshSubsets.h"

#include "AssemblerLib/SerialVectorMatrixBuilder.h"

template <typename Builder>
class AssemblerLibSerialVectorMatrixBuilder : public ::testing::Test
{
    public:
    typedef MeshLib::MeshItemType MeshItemType;
    typedef MeshLib::Location Location;
    typedef AssemblerLib::MeshComponentMap MeshComponentMap;

    typedef typename Builder::VectorType VectorType;
    typedef typename Builder::MatrixType MatrixType;

    public:
    AssemblerLibSerialVectorMatrixBuilder()
        : mesh(nullptr), nodesSubset(nullptr), cmap(nullptr)
    {
        mesh = MeshLib::MeshGenerator::generateLineMesh(1.0, mesh_size);
        nodesSubset = new MeshLib::MeshSubset(*mesh, mesh->getNodes());

        // Add two components both based on the same nodesSubset.
        components.emplace_back(new MeshLib::MeshSubsets(nodesSubset));
        components.emplace_back(new MeshLib::MeshSubsets(nodesSubset));

        cmap = new MeshComponentMap(components,
            AssemblerLib::ComponentOrder::BY_COMPONENT);
    }

    ~AssemblerLibSerialVectorMatrixBuilder()
    {
        delete cmap;
        std::remove_if(components.begin(), components.end(),
            [](MeshLib::MeshSubsets* p) { delete p; return true; });
        delete nodesSubset;
        delete mesh;
    }

    static std::size_t const mesh_size = 9;
    MeshLib::Mesh const* mesh;
    MeshLib::MeshSubset const* nodesSubset;

    std::vector<MeshLib::MeshSubsets*> components;
    MeshComponentMap const* cmap;
};

TYPED_TEST_CASE_P(AssemblerLibSerialVectorMatrixBuilder);

TYPED_TEST_P(AssemblerLibSerialVectorMatrixBuilder, createVector)
{
    typedef typename TestFixture::VectorType V;
    typedef TypeParam Builder;
    V* v = Builder::createVector(this->cmap->size());

    ASSERT_TRUE(v != nullptr);
    ASSERT_EQ(this->cmap->size(), v->size());

    delete v;
}

TYPED_TEST_P(AssemblerLibSerialVectorMatrixBuilder, createMatrix)
{
    typedef typename TestFixture::MatrixType M;
    typedef TypeParam Builder;
    M* m = Builder::createMatrix(this->cmap->size());

    ASSERT_TRUE(m != nullptr);
    ASSERT_EQ(this->cmap->size(), m->getNRows());
    ASSERT_EQ(this->cmap->size(), m->getNCols());

    delete m;
}

REGISTER_TYPED_TEST_CASE_P(AssemblerLibSerialVectorMatrixBuilder,
    createVector, createMatrix);

#include "MathLib/LinAlg/Dense/DenseVector.h"
#include "MathLib/LinAlg/Dense/GlobalDenseMatrix.h"

#ifdef USE_LIS
#include "MathLib/LinAlg/Lis/LisVector.h"
#include "MathLib/LinAlg/Lis/LisMatrix.h"
#endif  // USE_LIS

#ifdef USE_PETSC
#include "MathLib/LinAlg/PETSc/PETScVector.h"
#include "MathLib/LinAlg/PETSc/PETScMatrix.h"
#endif  // USE_PETSC

typedef ::testing::Types
    < AssemblerLib::SerialVectorMatrixBuilder<
        MathLib::GlobalDenseMatrix<double>, MathLib::DenseVector<double>>
#ifdef USE_LIS
    , AssemblerLib::SerialVectorMatrixBuilder<
        MathLib::LisMatrix, MathLib::LisVector>
#endif  // USE_LIS
#ifdef USE_PETSC
    , AssemblerLib::SerialVectorMatrixBuilder<
        MathLib::PETScMatrix, MathLib::PETScVector>
#endif  // USE_PETSC
    > TestTypes;

INSTANTIATE_TYPED_TEST_CASE_P(templated, AssemblerLibSerialVectorMatrixBuilder,
    TestTypes);

