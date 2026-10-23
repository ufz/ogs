// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <Eigen/Core>
#include <numeric>
#include <vector>

#include "MathLib/LinAlg/FinalizeMatrixAssembly.h"
#include "MathLib/LinAlg/FinalizeVectorAssembly.h"
#include "MathLib/LinAlg/LinAlg.h"
#include "MathLib/LinAlg/UnifiedMatrixSetters.h"

#if defined(USE_PETSC)
#include "MathLib/LinAlg/PETSc/PETScMatrix.h"
#include "MathLib/LinAlg/PETSc/PETScVector.h"
#include "Tests/MathLib/PETScSparsityPatternUtils.h"
#else
#include "MathLib/LinAlg/Eigen/EigenMatrix.h"
#include "MathLib/LinAlg/Eigen/EigenVector.h"
#endif

namespace
{
#if defined(USE_PETSC)
using TestMatrix = MathLib::PETScMatrix;
using TestVector = MathLib::PETScVector;

TestMatrix createTestMatrix(PetscInt const n)
{
    // Serial test: a single rank owns the whole dense n x n block.
    return TestMatrix(
        n, PETScSparsityPatternUtils::denseBlockSparsityPattern(n, n, n, 0));
}

std::vector<double> getVectorEntries(TestVector& v)
{
    MathLib::finalizeVectorAssembly(v);
    std::vector<double> values(v.size());
    v.getGlobalVector(values);
    return values;
}

Eigen::MatrixXd getMatrixEntries(TestMatrix& m)
{
    MathLib::finalizeMatrixAssembly(m);
    auto const rows = m.getNumberOfRows();
    auto const cols = m.getNumberOfColumns();

    std::vector<PetscInt> row_idcs(rows);
    std::vector<PetscInt> col_idcs(cols);
    std::iota(row_idcs.begin(), row_idcs.end(), 0);
    std::iota(col_idcs.begin(), col_idcs.end(), 0);

    // MatGetValues returns the values row-major.
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
        dense(rows, cols);
    MatGetValues(m.getRawMatrix(), rows, row_idcs.data(), cols, col_idcs.data(),
                 dense.data());
    return dense;
}
#else
using TestMatrix = MathLib::EigenMatrix;
using TestVector = MathLib::EigenVector;

TestMatrix createTestMatrix(TestMatrix::IndexType const n)
{
    return TestMatrix(n);
}

std::vector<double> getVectorEntries(TestVector& v)
{
    auto const& raw = v.getRawVector();
    return {raw.data(), raw.data() + raw.size()};
}

Eigen::MatrixXd getMatrixEntries(TestMatrix& m)
{
    return Eigen::MatrixXd(m.getRawMatrix());
}
#endif

TestVector createZeroTestVector(std::size_t const n)
{
    TestVector v(n);
    MathLib::LinAlg::set(v, 0.0);
    MathLib::finalizeVectorAssembly(v);
    return v;
}
}  // namespace

TEST(MathLibUnifiedMatrixSetters, SetVectorFromInitializerList)
{
    auto v = createZeroTestVector(4);

    MathLib::setVector(v, {1.0, -2.0, 3.5, 0.25});

    EXPECT_EQ((std::vector<double>{1.0, -2.0, 3.5, 0.25}), getVectorEntries(v));
}

TEST(MathLibUnifiedMatrixSetters, SetVectorFromInitializerListOverwrites)
{
    auto v = createZeroTestVector(3);
    MathLib::setVector(v, {1.0, 2.0, 3.0});
    MathLib::finalizeVectorAssembly(v);

    MathLib::setVector(v, {4.0, 5.0, 6.0});

    EXPECT_EQ((std::vector<double>{4.0, 5.0, 6.0}), getVectorEntries(v));
}

TEST(MathLibUnifiedMatrixSetters, SetVectorSingleEntry)
{
    auto v = createZeroTestVector(4);
    MathLib::setVector(v, {1.0, 2.0, 3.0, 4.0});
    MathLib::finalizeVectorAssembly(v);

    MathLib::setVector(v, 2, -7.0);

    // Only the addressed entry changes; it is overwritten, not added to.
    EXPECT_EQ((std::vector<double>{1.0, 2.0, -7.0, 4.0}), getVectorEntries(v));
}

TEST(MathLibUnifiedMatrixSetters, SetVectorMultipleEntries)
{
    auto v = createZeroTestVector(5);
    MathLib::setVector(v, {1.0, 2.0, 3.0, 4.0, 5.0});
    MathLib::finalizeVectorAssembly(v);

    // Unordered indices: each value must go to its paired index.
    MathLib::setVector(v, {4, 0, 2}, {-4.0, -0.5, -2.0});

    EXPECT_EQ((std::vector<double>{-0.5, 2.0, -2.0, 4.0, -4.0}),
              getVectorEntries(v));
}

TEST(MathLibUnifiedMatrixSetters, SetVectorMultipleEntriesFromStdVectors)
{
    auto v = createZeroTestVector(4);

    std::vector<MathLib::MatrixVectorTraits<TestVector>::Index> const indices{
        3, 1};
    std::vector<double> const values{8.0, -1.5};
    MathLib::setVector(v, indices, values);

    EXPECT_EQ((std::vector<double>{0.0, -1.5, 0.0, 8.0}), getVectorEntries(v));
}

TEST(MathLibUnifiedMatrixSetters, SetMatrixFromInitializerList)
{
    auto m = createTestMatrix(2);

    // Values are given row by row.
    MathLib::setMatrix(m, {1.0, 2.0, 3.0, 4.0});

    Eigen::Matrix2d expected;
    expected << 1.0, 2.0, 3.0, 4.0;
    EXPECT_EQ(expected, getMatrixEntries(m));
}

TEST(MathLibUnifiedMatrixSetters, SetMatrixFromInitializerListOverwrites)
{
    auto m = createTestMatrix(2);
    MathLib::setMatrix(m, {1.0, 2.0, 3.0, 4.0});
    MathLib::finalizeMatrixAssembly(m);

    // Previously non-zero entries must be reset, not kept or accumulated.
    MathLib::setMatrix(m, {0.0, -1.0, 5.0, 0.0});

    Eigen::Matrix2d expected;
    expected << 0.0, -1.0, 5.0, 0.0;
    EXPECT_EQ(expected, getMatrixEntries(m));
}

TEST(MathLibUnifiedMatrixSetters, SetMatrixFromEigenMatrix)
{
    auto m = createTestMatrix(3);
    MathLib::setMatrix(m, {9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0});
    MathLib::finalizeMatrixAssembly(m);

    // Asymmetric, so a row-/column-major mix-up would be detected.
    Eigen::MatrixXd values(3, 3);
    values << 1.0, 2.0, 0.0,  //
        0.0, 4.0, 5.0,        //
        7.0, 0.0, 9.0;

    MathLib::setMatrix(m, values);

    EXPECT_EQ(values, getMatrixEntries(m));
}

TEST(MathLibUnifiedMatrixSetters, AddToMatrix)
{
    auto m = createTestMatrix(2);
    MathLib::setMatrix(m, {1.0, 2.0, 3.0, 4.0});
    MathLib::finalizeMatrixAssembly(m);

    MathLib::addToMatrix(m, {10.0, 0.0, -3.0, 0.5});

    Eigen::Matrix2d expected;
    expected << 11.0, 2.0, 0.0, 4.5;
    EXPECT_EQ(expected, getMatrixEntries(m));
}
