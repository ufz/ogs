/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <gtest/gtest.h>
#include <limits>
#include <logog/include/logog.hpp>

#include "MathLib/LinAlg/LinAlg.h"
#include "MathLib/LinAlg/MatrixSpecifications.h"
#include "MathLib/LinAlg/MatrixVectorTraits.h"
#include "MeshLib/MeshGenerators/MeshGenerator.h"
#include "NumLib/DOF/LocalToGlobalIndexMap.h"
#include "NumLib/NumericsConfig.h"
#include "NumLib/DOF/DOFTableUtil.h"
#include "Tests/VectorUtils.h"

static std::vector<std::unique_ptr<MeshLib::MeshSubsets>> createMeshSubsets(
    MeshLib::MeshSubset const& mesh_subset, unsigned const num_components)
{
    std::vector<std::unique_ptr<MeshLib::MeshSubsets>> mesh_subsets;
    for (unsigned i=0; i<num_components; ++i)
        mesh_subsets.emplace_back(
            new MeshLib::MeshSubsets{&mesh_subset});

    return mesh_subsets;
}

struct DOFTableData
{
    DOFTableData(unsigned const num_components,
                 NumLib::ComponentOrder const order)
        : mesh(MeshLib::MeshGenerator::generateRegularHexMesh(
              mesh_length, mesh_elements_in_each_direction)),
          mesh_subset_all_nodes(*mesh, &mesh->getNodes()),
          dof_table(createMeshSubsets(mesh_subset_all_nodes, num_components),
                    order)
    {
    }

    std::unique_ptr<MeshLib::Mesh> const mesh;
    MeshLib::MeshSubset const mesh_subset_all_nodes;
    NumLib::LocalToGlobalIndexMap const dof_table;

private:
    static const double mesh_length;
    static const std::size_t mesh_elements_in_each_direction;
};
const double DOFTableData::mesh_length = 1;
const std::size_t DOFTableData::mesh_elements_in_each_direction = 5;

template <typename AccumulateCallback, typename AccumulateFinishCallback>
void do_test(unsigned const num_components,
             MathLib::VecNormType const norm_type,
             double const tolerance,
             AccumulateCallback&& accumulate_cb,
             AccumulateFinishCallback&& accumulate_finish_cb)
{
    /* The idea behind this test is that the accumulated value of all component
     * norms must be equal to the 1-/2-/infinity-norm of the whole global
     * vector.
     */
    using CO = NumLib::ComponentOrder;
    for (auto order : {CO::BY_COMPONENT, CO::BY_LOCATION})
    {
        DOFTableData dtd(num_components, order);

        MathLib::MatrixSpecifications mat_specs(
            dtd.dof_table.dofSizeWithoutGhosts(),
            dtd.dof_table.dofSizeWithoutGhosts(),
            &dtd.dof_table.getGhostIndices(),
            nullptr);

        auto x =
            MathLib::MatrixVectorTraits<GlobalVector>::newInstance(mat_specs);

        fillVectorRandomly(*x);

        auto const total_norm = MathLib::LinAlg::norm(*x, norm_type);

        double compwise_total_norm = 0.0;
        for (unsigned comp = 0; comp < num_components; ++comp) {
            compwise_total_norm = accumulate_cb(
                compwise_total_norm,
                NumLib::norm(*x, comp, norm_type, dtd.dof_table, *dtd.mesh));
        }
        compwise_total_norm = accumulate_finish_cb(compwise_total_norm);

        EXPECT_NEAR(total_norm, compwise_total_norm, tolerance);
    }
}

#ifndef USE_PETSC
TEST(NumLib, ComponentNormSingleComponent)
#else
TEST(NumLib, DISABLED_ComponentNormSingleComponent)
#endif
{
    unsigned const num_components = 1;
    auto const tolerance = 800 * std::numeric_limits<double>::epsilon();

    using VNT = MathLib::VecNormType;
    for (auto norm_type : {VNT::NORM1, VNT::NORM2, VNT::INFINITY_N}) {
        DBUG("norm type: %s.",
             MathLib::convertVecNormTypeToString(norm_type).c_str());
        do_test(num_components, norm_type, tolerance,
                [](double /*n_total*/, double n) { return n; },
                [](double n_total) { return n_total; });
    }
}

#ifndef USE_PETSC
TEST(NumLib, ComponentNormMultiComponent1)
#else
TEST(NumLib, DISABLED_ComponentNormMultiComponent1)
#endif
{
    unsigned const num_components = 3;
    auto const norm_type = MathLib::VecNormType::NORM1;
    auto const tolerance =
        num_components * 600 * std::numeric_limits<double>::epsilon();

    do_test(num_components, norm_type, tolerance,
            [](double n_total, double n) { return n_total + n; },
            [](double n_total) { return n_total; });
}

#ifndef USE_PETSC
TEST(NumLib, ComponentNormMultiComponent2)
#else
TEST(NumLib, DISABLED_ComponentNormMultiComponent2)
#endif
{
    unsigned const num_components = 3;
    auto const norm_type = MathLib::VecNormType::NORM2;
    auto const tolerance =
        num_components * 20 * std::numeric_limits<double>::epsilon();

    do_test(num_components, norm_type, tolerance,
            [](double n_total, double n) { return n_total + n*n; },
            [](double n_total) { return std::sqrt(n_total); });
}

#ifndef USE_PETSC
TEST(NumLib, ComponentNormMultiComponentInfinity)
#else
TEST(NumLib, DISABLED_ComponentNormMultiComponentInfinity)
#endif
{
    unsigned const num_components = 3;
    auto const norm_type = MathLib::VecNormType::INFINITY_N;
    auto const tolerance =
        num_components * 1 * std::numeric_limits<double>::epsilon();

    do_test(num_components, norm_type, tolerance,
            [](double n_total, double n) { return std::max(n_total, n); },
            [](double n_total) { return n_total; });
}
