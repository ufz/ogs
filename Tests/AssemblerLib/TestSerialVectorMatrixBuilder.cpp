/**
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <gtest/gtest.h>

#include "AssemblerLib/MeshComponentMap.h"

#include "MeshLib/MeshGenerators/MeshGenerator.h"
#include "MeshLib/MeshSubsets.h"

#include "AssemblerLib/SerialDenseVectorMatrixBuilder.h"

#ifdef USE_LIS
#include "AssemblerLib/SerialLisVectorMatrixBuilder.h"
#endif  // USE_LIS

template <typename SerialBuilder>
class AssemblerLibSerialVectorMatrixBuilder : public ::testing::Test
{
    public:
    typedef MeshLib::MeshItemType MeshItemType;
    typedef MeshLib::Location Location;
    typedef AssemblerLib::MeshComponentMap MeshComponentMap;

    typedef typename SerialBuilder::VectorType VectorType;
    typedef typename SerialBuilder::MatrixType MatrixType;

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
    typedef TypeParam SerialBuilder;
    V* v = SerialBuilder::createVector(this->cmap->size());

    ASSERT_TRUE(v != nullptr);
    ASSERT_EQ(this->cmap->size(), v->size());

    delete v;
}

TYPED_TEST_P(AssemblerLibSerialVectorMatrixBuilder, createMatrix)
{
    typedef typename TestFixture::MatrixType M;
    typedef TypeParam SerialBuilder;
    M* v = SerialBuilder::createMatrix(this->cmap->size());

    ASSERT_TRUE(v != nullptr);
    ASSERT_EQ(this->cmap->size(), v->getNRows());
    ASSERT_EQ(this->cmap->size(), v->getNCols());

    delete v;
}

REGISTER_TYPED_TEST_CASE_P(AssemblerLibSerialVectorMatrixBuilder,
    createVector, createMatrix);


typedef ::testing::Types
    < AssemblerLib::SerialDenseVectorMatrixBuilder
#ifdef USE_LIS
    , AssemblerLib::SerialLisVectorMatrixBuilder
#endif  // USE_LIS
    > TestTypes;

INSTANTIATE_TYPED_TEST_CASE_P(templated, AssemblerLibSerialVectorMatrixBuilder,
    TestTypes);

