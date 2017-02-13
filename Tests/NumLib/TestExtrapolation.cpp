/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <gtest/gtest.h>

#include "MathLib/LinAlg/LinAlg.h"

#include "MeshLib/MeshGenerators/MeshGenerator.h"
#include "MeshLib/IO/writeMeshToFile.h"

#include "NumLib/DOF/DOFTableUtil.h"
#include "NumLib/DOF/MatrixProviderUser.h"
#include "NumLib/Extrapolation/ExtrapolatableElementCollection.h"
#include "NumLib/Extrapolation/Extrapolator.h"
#include "NumLib/Extrapolation/LocalLinearLeastSquaresExtrapolator.h"
#include "NumLib/Fem/FiniteElement/TemplateIsoparametric.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "NumLib/Function/Interpolation.h"
#include "NumLib/NumericsConfig.h"

#include "ProcessLib/Utils/LocalDataInitializer.h"
#include "ProcessLib/Utils/CreateLocalAssemblers.h"
#include "ProcessLib/Utils/InitShapeMatrices.h"

#include "Tests/VectorUtils.h"

namespace
{
template <typename ShapeMatrices>
void interpolateNodalValuesToIntegrationPoints(
    std::vector<double> const& local_nodal_values,
    std::vector<ShapeMatrices, Eigen::aligned_allocator<ShapeMatrices>> const&
        shape_matrices,
    std::vector<double>& interpolated_values)
{
    for (unsigned ip=0; ip<shape_matrices.size(); ++ip)
    {
        NumLib::shapeFunctionInterpolate(
            local_nodal_values, shape_matrices[ip].N, interpolated_values[ip]);
    }
}

} // anonymous namespace

class LocalAssemblerDataInterface : public NumLib::ExtrapolatableElement
{
public:
    virtual void interpolateNodalValuesToIntegrationPoints(
        std::vector<double> const& local_nodal_values) = 0;

    virtual std::vector<double> const& getStoredQuantity(
        std::vector<double>& /*cache*/) const = 0;

    virtual std::vector<double> const& getDerivedQuantity(
        std::vector<double>& cache) const = 0;
};

using IntegrationPointValuesMethod = std::vector<double> const& (
    LocalAssemblerDataInterface::*)(std::vector<double>&)const;

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
class LocalAssemblerData : public LocalAssemblerDataInterface
{
    using ShapeMatricesType = ShapeMatrixPolicyType<ShapeFunction, GlobalDim>;
    using ShapeMatrices = typename ShapeMatricesType::ShapeMatrices;

public:
    LocalAssemblerData(MeshLib::Element const& e,
                       std::size_t const /*local_matrix_size*/,
                       bool is_axially_symmetric,
                       unsigned const integration_order)
        : _shape_matrices(
              ProcessLib::initShapeMatrices<ShapeFunction, ShapeMatricesType,
                                            IntegrationMethod, GlobalDim>(
                  e, is_axially_symmetric,
                  IntegrationMethod{integration_order})),
          _int_pt_values(_shape_matrices.size())
    {
    }

    Eigen::Map<const Eigen::RowVectorXd>
    getShapeMatrix(const unsigned integration_point) const override
    {
        auto const& N = _shape_matrices[integration_point].N;

        // assumes N is stored contiguously in memory
        return Eigen::Map<const Eigen::RowVectorXd>(N.data(), N.size());
    }

    std::vector<double> const& getStoredQuantity(
        std::vector<double>& /*cache*/) const override
    {
        return _int_pt_values;
    }

    std::vector<double> const& getDerivedQuantity(
        std::vector<double>& cache) const override
    {
        cache.clear();
        for (auto value : _int_pt_values)
            cache.push_back(2.0 * value);
        return cache;
    }

    void interpolateNodalValuesToIntegrationPoints(
        std::vector<double> const& local_nodal_values) override
    {
        ::interpolateNodalValuesToIntegrationPoints(
            local_nodal_values, _shape_matrices, _int_pt_values);
    }

private:
    std::vector<ShapeMatrices, Eigen::aligned_allocator<ShapeMatrices>>
        _shape_matrices;
    std::vector<double> _int_pt_values;
};

class TestProcess
{
public:
    using LocalAssembler = LocalAssemblerDataInterface;

    using ExtrapolatorInterface = NumLib::Extrapolator;
    using ExtrapolatorImplementation =
        NumLib::LocalLinearLeastSquaresExtrapolator;

    TestProcess(MeshLib::Mesh const& mesh, unsigned const integration_order)
        : _integration_order(integration_order)
        , _mesh_subset_all_nodes(mesh, &mesh.getNodes())
    {
        std::vector<std::unique_ptr<MeshLib::MeshSubsets>> all_mesh_subsets;
        all_mesh_subsets.emplace_back(
            new MeshLib::MeshSubsets{&_mesh_subset_all_nodes});

        _dof_table.reset(new NumLib::LocalToGlobalIndexMap(
            std::move(all_mesh_subsets), NumLib::ComponentOrder::BY_COMPONENT));

        // Passing _dof_table works, because this process has only one variable
        // and the variable has exactly one component.
        _extrapolator.reset(new ExtrapolatorImplementation(*_dof_table));

        // createAssemblers(mesh);
        ProcessLib::createLocalAssemblers<LocalAssemblerData>(
            mesh.getDimension(), mesh.getElements(), *_dof_table, 1,
            _local_assemblers, mesh.isAxiallySymmetric(), _integration_order);
    }

    void interpolateNodalValuesToIntegrationPoints(
        GlobalVector const& global_nodal_values) const
    {
        auto cb = [](std::size_t id, LocalAssembler& loc_asm,
                     NumLib::LocalToGlobalIndexMap const& dof_table,
                     GlobalVector const& x) {
            auto const indices = NumLib::getIndices(id, dof_table);
            auto const local_x = x.get(indices);

            loc_asm.interpolateNodalValuesToIntegrationPoints(local_x);
        };

        GlobalExecutor::executeDereferenced(cb, _local_assemblers, *_dof_table,
                                            global_nodal_values);
    }

    std::pair<GlobalVector const*, GlobalVector const*> extrapolate(
        IntegrationPointValuesMethod method) const
    {
        auto const extrapolatables =
            NumLib::makeExtrapolatable(_local_assemblers, method);

        _extrapolator->extrapolate(extrapolatables);
        _extrapolator->calculateResiduals(extrapolatables);

        return {&_extrapolator->getNodalValues(),
                &_extrapolator->getElementResiduals()};
    }

private:
    unsigned const _integration_order;

    MeshLib::MeshSubset _mesh_subset_all_nodes;
    std::unique_ptr<NumLib::LocalToGlobalIndexMap> _dof_table;

    std::vector<std::unique_ptr<LocalAssembler>> _local_assemblers;

    std::unique_ptr<ExtrapolatorInterface> _extrapolator;
};

void extrapolate(TestProcess const& pcs, IntegrationPointValuesMethod method,
                 GlobalVector const& expected_extrapolated_global_nodal_values,
                 std::size_t const nnodes, std::size_t const nelements)
{
    namespace LinAlg = MathLib::LinAlg;

    auto const tolerance_dx  = 30.0 * std::numeric_limits<double>::epsilon();
    auto const tolerance_res = 15.0 * std::numeric_limits<double>::epsilon();

    auto const result = pcs.extrapolate(method);
    auto const& x_extra = *result.first;
    auto const& residual = *result.second;

    ASSERT_EQ(nnodes,    x_extra.size());
    ASSERT_EQ(nelements, residual.size());

    auto const res_norm = LinAlg::normMax(residual);
    DBUG("maximum norm of residual: %g", res_norm);
    EXPECT_GT(tolerance_res, res_norm);

    auto delta_x = MathLib::MatrixVectorTraits<GlobalVector>::newInstance(
                expected_extrapolated_global_nodal_values);
    LinAlg::axpy(*delta_x, -1.0, x_extra); // delta_x = x_expected - x_extra

    auto const dx_norm = LinAlg::normMax(*delta_x);
    DBUG("maximum norm of delta x:  %g", dx_norm);
    EXPECT_GT(tolerance_dx, dx_norm);
}

#ifndef USE_PETSC
TEST(NumLib, Extrapolation)
#else
TEST(NumLib, DISABLED_Extrapolation)
#endif
{
    /* In this test a random vector x of nodal values is created.
     * This vector is interpolated to the integration points using each
     * element's the shape functions.
     * Afterwards the integration point values y are extrapolated to the mesh
     * nodes again.
     * Since y have been obtained from x via interpolation, it is expected, that
     * the interpolation result nearly exactly matches the original nodal values
     * x.
     */

    const double mesh_length = 1.0;
    const std::size_t mesh_elements_in_each_direction = 5;

    // generate mesh
    std::unique_ptr<MeshLib::Mesh> mesh(
                MeshLib::MeshGenerator::generateRegularHexMesh(
                    mesh_length, mesh_elements_in_each_direction));

    for (unsigned integration_order : {2, 3, 4})
    {
        namespace LinAlg = MathLib::LinAlg;

        auto const nnodes    = mesh->getNumberOfNodes();
        auto const nelements = mesh->getNumberOfElements();
        DBUG("number of nodes: %lu, number of elements: %lu", nnodes, nelements);

        TestProcess pcs(*mesh, integration_order);

        // generate random nodal values
        MathLib::MatrixSpecifications spec{nnodes, nnodes, nullptr, nullptr};
        auto x = MathLib::MatrixVectorTraits<GlobalVector>::newInstance(spec);

        fillVectorRandomly(*x);

        pcs.interpolateNodalValuesToIntegrationPoints(*x);

        // test extrapolation of a quantity that is stored in the local
        // assembler
        extrapolate(pcs, &LocalAssemblerDataInterface::getStoredQuantity, *x,
                    nnodes, nelements);

        // expect 2*x as extraplation result for derived quantity
        auto two_x = MathLib::MatrixVectorTraits<GlobalVector>::newInstance(*x);
        LinAlg::axpy(*two_x, 1.0, *x); // two_x = x + x

        // test extrapolation of a quantity that is derived from some
        // integration point values
        extrapolate(pcs, &LocalAssemblerDataInterface::getDerivedQuantity,
                    *two_x, nnodes, nelements);
    }
}
