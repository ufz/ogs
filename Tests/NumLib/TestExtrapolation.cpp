/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <random>
#include <gtest/gtest.h>

#include "NumLib/DOF/MatrixProviderUser.h"
#include "MathLib/LinAlg/MatrixVectorTraits.h"
#include "MathLib/LinAlg/UnifiedMatrixSetters.h"
#include "NumLib/Assembler/VectorMatrixAssembler.h"

#include "MathLib/LinAlg/BLAS.h"

#include "MeshLib/MeshGenerators/MeshGenerator.h"

#include "NumLib/Extrapolation/Extrapolator.h"
#include "NumLib/Extrapolation/LocalLinearLeastSquaresExtrapolator.h"
#include "NumLib/Fem/FiniteElement/TemplateIsoparametric.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"

#include "NumLib/NumericsConfig.h"
#include "ProcessLib/Utils/LocalDataInitializer.h"
#include "ProcessLib/Utils/CreateLocalAssemblers.h"
#include "ProcessLib/Utils/InitShapeMatrices.h"


namespace
{

template<typename ShapeMatrices>
void interpolateNodalValuesToIntegrationPoints(
        std::vector<double> const& local_nodal_values,
        std::vector<ShapeMatrices, Eigen::aligned_allocator<ShapeMatrices>> const& shape_matrices,
        std::vector<double>& interpolated_values)
{
    for (unsigned ip=0; ip<shape_matrices.size(); ++ip)
    {
        NumLib::shapeFunctionInterpolate(
            local_nodal_values, shape_matrices[ip].N, interpolated_values[ip]);
    }
}

template<typename Vector>
void fillVectorRandomly(Vector& x)
{
    std::random_device rd;
    std::mt19937 random_number_generator(rd());
    std::uniform_real_distribution<double> rnd;
    std::vector<GeoLib::Point> pnts;

    using Index = typename MathLib::MatrixVectorTraits<Vector>::Index;
    Index const size = x.size();

    for (Index i=0; i<size; ++i) {
        MathLib::setVector(x, i, rnd(random_number_generator));
    }
}

}

enum class IntegrationPointValue
{
    StoredQuantity, // some quantity acutally stored in the local assembler
    DerivedQuantity // a quantity computed for each integration point on-the-fly
};

template<typename GlobalMatrix, typename GlobalVector>
class LocalAssemblerDataInterface
        : public NumLib::Extrapolatable<GlobalVector, IntegrationPointValue>
{
public:
    virtual void interpolateNodalValuesToIntegrationPoints(
            std::vector<double> const& local_nodal_values,
            IntegrationPointValue const property) = 0;
};

template<typename ShapeFunction,
         typename IntegrationMethod,
         typename GlobalMatrix,
         typename GlobalVector,
         unsigned GlobalDim>
class LocalAssemblerData
        : public LocalAssemblerDataInterface<GlobalMatrix, GlobalVector>
{
    using ShapeMatricesType = ShapeMatrixPolicyType<ShapeFunction, GlobalDim>;
    using ShapeMatrices = typename ShapeMatricesType::ShapeMatrices;

public:
    LocalAssemblerData(MeshLib::Element const& e,
                       std::size_t const /*local_matrix_size*/,
                       unsigned const integration_order)
        : _shape_matrices(
              ProcessLib::initShapeMatrices<ShapeFunction, ShapeMatricesType,
                  IntegrationMethod, GlobalDim>(e, integration_order))
        , _int_pt_values(_shape_matrices.size())
    {}

    Eigen::Map<const Eigen::VectorXd>
    getShapeMatrix(const unsigned integration_point) const override
    {
        auto const& N = _shape_matrices[integration_point].N;

        // assumes N is stored contiguously in memory
        return Eigen::Map<const Eigen::VectorXd>(N.data(), N.size());
    }

    std::vector<double> const&
    getIntegrationPointValues(IntegrationPointValue const property,
                              std::vector<double>& cache) const override
    {
        switch (property)
        {
        case IntegrationPointValue::StoredQuantity:
            return _int_pt_values;
        case IntegrationPointValue::DerivedQuantity:
            cache.clear();
            for (auto value : _int_pt_values)
                cache.push_back(2.0 * value);
            return cache;
        }

        OGS_FATAL("");
    }

    void interpolateNodalValuesToIntegrationPoints(
            std::vector<double> const& local_nodal_values,
            IntegrationPointValue const property) override
    {
        switch (property)
        {
        case IntegrationPointValue::StoredQuantity:
            ::interpolateNodalValuesToIntegrationPoints(
                        local_nodal_values, _shape_matrices, _int_pt_values);
            break;
        case IntegrationPointValue::DerivedQuantity:
            break;
        }
    }

private:
    std::vector<ShapeMatrices, Eigen::aligned_allocator<ShapeMatrices>> _shape_matrices;
    std::vector<double>        _int_pt_values;
};

template<typename GlobalSetup>
class TestProcess
{
public:
    using GlobalMatrix = typename GlobalSetup::MatrixType;
    using GlobalVector = typename GlobalSetup::VectorType;

    using LocalAssembler = LocalAssemblerDataInterface<GlobalMatrix, GlobalVector>;
    using GlobalAssembler = NumLib::VectorMatrixAssembler<
        GlobalMatrix, GlobalVector, LocalAssembler,
        // The exact tag does not matter here.
        NumLib::ODESystemTag::FirstOrderImplicitQuasilinear>;

    using ExtrapolatorInterface =
        NumLib::Extrapolator<GlobalVector, IntegrationPointValue, LocalAssembler>;
    using ExtrapolatorImplementation =
        NumLib::LocalLinearLeastSquaresExtrapolator<
            GlobalVector, IntegrationPointValue, LocalAssembler>;

    TestProcess(MeshLib::Mesh const& mesh, unsigned const integration_order)
        : _integration_order(integration_order)
        , _mesh_subset_all_nodes(mesh, &mesh.getNodes())
    {
        std::vector<std::unique_ptr<MeshLib::MeshSubsets>> all_mesh_subsets;
        all_mesh_subsets.emplace_back(
                    new MeshLib::MeshSubsets{&_mesh_subset_all_nodes});

        _dof_table.reset(new NumLib::LocalToGlobalIndexMap(
              std::move(all_mesh_subsets),
              NumLib::ComponentOrder::BY_COMPONENT));

        // Passing _dof_table works, because this process has only one variable
        // and the variable has exactly one component.
        _extrapolator.reset(new ExtrapolatorImplementation(
            MathLib::MatrixSpecifications{_dof_table->dofSizeWithoutGhosts(),
                                          _dof_table->dofSizeWithoutGhosts(),
                                          &_dof_table->getGhostIndices(),
                                          nullptr},
            *_dof_table));

        createAssemblers(mesh);
    }

    void interpolateNodalValuesToIntegrationPoints(
            GlobalVector const& global_nodal_values,
            IntegrationPointValue const property) const
    {
        auto cb = [this, property](
            std::size_t id, LocalAssembler& loc_asm, GlobalVector const& x)
        {
            auto inner_cb = [&loc_asm, property](
                std::vector<double> const& local_x,
                NumLib::LocalToGlobalIndexMap::RowColumnIndices const&
            ) {
                loc_asm.interpolateNodalValuesToIntegrationPoints(local_x, property);
            };

            _global_assembler->passLocalVector(inner_cb, id, x);
        };

        GlobalSetup::executeDereferenced(
                    cb, _local_assemblers, global_nodal_values);
    }

    std::pair<GlobalVector const*, GlobalVector const*>
    extrapolate(IntegrationPointValue const property) const
    {
        _extrapolator->extrapolate(_local_assemblers, property);
        _extrapolator->calculateResiduals(_local_assemblers, property);

        return { &_extrapolator->getNodalValues(),
                 &_extrapolator->getElementResiduals() };
    }

private:
    void createAssemblers(MeshLib::Mesh const& mesh)
    {
        _global_assembler.reset(new GlobalAssembler(*_dof_table));
        switch (mesh.getDimension())
        {
        case 1:  createLocalAssemblers<1>(mesh); break;
        case 2:  createLocalAssemblers<2>(mesh); break;
        case 3:  createLocalAssemblers<3>(mesh); break;
        default: assert(false);
        }
    }

    template<unsigned GlobalDim>
    void createLocalAssemblers(MeshLib::Mesh const& mesh)
    {
        using LocalDataInitializer = ProcessLib::LocalDataInitializer<
            LocalAssembler, LocalAssemblerData,
            GlobalMatrix, GlobalVector, GlobalDim>;

        _local_assemblers.resize(mesh.getNumberOfElements());

        LocalDataInitializer initializer(*_dof_table);

        DBUG("Calling local assembler builder for all mesh elements.");
        GlobalSetup::transformDereferenced(
                initializer,
                mesh.getElements(),
                _local_assemblers,
                _integration_order);
    }

    unsigned const _integration_order;

    MeshLib::MeshSubset _mesh_subset_all_nodes;
    std::unique_ptr<NumLib::LocalToGlobalIndexMap> _dof_table;

    std::unique_ptr<GlobalAssembler> _global_assembler;
    std::vector<std::unique_ptr<LocalAssembler>> _local_assemblers;

    std::unique_ptr<ExtrapolatorInterface> _extrapolator;
};


template<typename GlobalSetup>
void extrapolate(TestProcess<GlobalSetup> const& pcs,
                 IntegrationPointValue property,
                 typename GlobalSetup::VectorType const&
                 expected_extrapolated_global_nodal_values,
                 std::size_t const nnodes, std::size_t const nelements)
{
    namespace BLAS = MathLib::BLAS;
    using GlobalVector = typename GlobalSetup::VectorType;

    auto const tolerance_dx  = 20.0 * std::numeric_limits<double>::epsilon();
    auto const tolerance_res =  4.0 * std::numeric_limits<double>::epsilon();

    auto const  result   = pcs.extrapolate(property);
    auto const& x_extra  = *result.first;
    auto const& residual = *result.second;

    ASSERT_EQ(nnodes,    x_extra.size());
    ASSERT_EQ(nelements, residual.size());

    auto const res_norm = BLAS::normMax(residual);
    DBUG("maximum norm of residual: %g", res_norm);
    EXPECT_GT(tolerance_res, res_norm);

    auto delta_x = MathLib::MatrixVectorTraits<GlobalVector>::newInstance(
                expected_extrapolated_global_nodal_values);
    BLAS::axpy(*delta_x, -1.0, x_extra); // delta_x = x_expected - x_extra

    auto const dx_norm = BLAS::normMax(*delta_x);
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

    for (unsigned integration_order : {2, 3, 4})
    {

        namespace BLAS = MathLib::BLAS;
        using GlobalSetup = GlobalSetupType;
        using GlobalVector = GlobalSetup::VectorType;

        const double mesh_length = 1.0;
        const double mesh_elements_in_each_direction = 5.0;

        // generate mesh
        std::unique_ptr<MeshLib::Mesh> mesh(
                    MeshLib::MeshGenerator::generateRegularHexMesh(
                        mesh_length, mesh_elements_in_each_direction));

        auto const nnodes    = mesh->getNumberOfNodes();
        auto const nelements = mesh->getNumberOfElements();
        DBUG("number of nodes: %lu, number of elements: %lu", nnodes, nelements);

        TestProcess<GlobalSetup> pcs(*mesh, integration_order);

        // generate random nodal values
        MathLib::MatrixSpecifications spec{nnodes, nnodes, nullptr, nullptr};
        auto x = MathLib::MatrixVectorTraits<GlobalVector>::newInstance(spec);

        fillVectorRandomly(*x);

        pcs.interpolateNodalValuesToIntegrationPoints(
                    *x, IntegrationPointValue::StoredQuantity);

        // test extrapolation of a quantity that is stored in the local assembler
        extrapolate(pcs, IntegrationPointValue::StoredQuantity,
                    *x, nnodes, nelements);

        // expect 2*x as extraplation result for derived quantity
        auto two_x = MathLib::MatrixVectorTraits<GlobalVector>::newInstance(*x);
        BLAS::axpy(*two_x, 1.0, *x); // two_x = x + x

        // test extrapolation of a quantity that is derived from some integration
        // point values
        extrapolate(pcs, IntegrationPointValue::DerivedQuantity,
                    *two_x, nnodes, nelements);
    }
}
