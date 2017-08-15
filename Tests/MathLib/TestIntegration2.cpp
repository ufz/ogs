/**
 * \author Norihiro Watanabe
 * \date   2013-08-29
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include <gtest/gtest.h>

#include <cmath>
#include <limits>

#include "BaseLib/BuildInfo.h"

#include "MathLib/Integration/GaussLegendreTet.h"

#include "MeshLib/Elements/Element.h"
#include "MeshLib/IO/VtkIO/VtuInterface.h"
#include "MeshLib/MeshSubsets.h"

#include "NumLib/DOF/LocalToGlobalIndexMap.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"

#include "ProcessLib/Utils/CreateLocalAssemblers.h"
#include "ProcessLib/Utils/InitShapeMatrices.h"

#include "Tests/VectorUtils.h"

template <typename Shp>
static std::array<double, 3> interpolateNodeCoordinates(
    MeshLib::Element const& e, Shp const& N)
{
    std::array<double, 3> res;

    auto const* const nodes = e.getNodes();
    auto node_coords = N;

    for (std::size_t d = 0; d < 3; ++d)
    {
        for (unsigned ip = 0; ip < N.size(); ++ip)
        {
            node_coords[ip] = (*nodes[ip])[d];
        }

        res[d] = N.dot(node_coords);
    }

    return res;
}

class LocalAssemblerDataInterface
{
public:
    using Function = std::function<double(std::array<double, 3> const&)>;

    virtual double integrate(Function const& f,
                             unsigned const integration_order) const = 0;

    virtual ~LocalAssemblerDataInterface() = default;
};

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
                       unsigned const /*integration_order*/)
        : _e(e)
    {
        if (is_axially_symmetric)
            OGS_FATAL("Only testing Cartesian meshes!");
    }

    double integrate(const Function& f,
                     unsigned const integration_order) const override
    {
        double integral = 0;

        auto const sms =
            ProcessLib::initShapeMatrices<ShapeFunction, ShapeMatricesType,
                                          IntegrationMethod, GlobalDim>(
                _e, false /*is_axially_symmetric*/,
                IntegrationMethod{integration_order});
        IntegrationMethod integration_method{integration_order};

        for (unsigned ip = 0; ip < sms.size(); ++ip)
        {
            auto const& N = sms[ip].N;
            auto const coords = interpolateNodeCoordinates(_e, N);
            auto const function_value = f(coords);
            integral += function_value * sms[ip].detJ *
                        integration_method.getWeightedPoint(ip).getWeight();
        }

        return integral;
    }

private:
    MeshLib::Element const& _e;
};

class TestProcess
{
public:
    using LocalAssembler = LocalAssemblerDataInterface;

    TestProcess(MeshLib::Mesh const& mesh, unsigned const integration_order)
        : _integration_order(integration_order),
          _mesh_subset_all_nodes(mesh, &mesh.getNodes())
    {
        std::vector<MeshLib::MeshSubsets> all_mesh_subsets;
        all_mesh_subsets.emplace_back(&_mesh_subset_all_nodes);

        _dof_table = std::make_unique<NumLib::LocalToGlobalIndexMap>(
            std::move(all_mesh_subsets), NumLib::ComponentOrder::BY_COMPONENT);

        // createAssemblers(mesh);
        ProcessLib::createLocalAssemblers<LocalAssemblerData>(
            mesh.getDimension(), mesh.getElements(), *_dof_table, 1,
            _local_assemblers, mesh.isAxiallySymmetric(), _integration_order);
    }

    double integrate(LocalAssembler::Function const& f) const
    {
        double integral = 0;
        for (auto const& la : _local_assemblers)
        {
            integral += la->integrate(f, _integration_order);
        }

        return integral;
    }

private:
    unsigned const _integration_order;

    MeshLib::MeshSubset _mesh_subset_all_nodes;
    std::unique_ptr<NumLib::LocalToGlobalIndexMap> _dof_table;

    std::vector<std::unique_ptr<LocalAssembler>> _local_assemblers;
};

LocalAssemblerDataInterface::Function getFConst(
    std::vector<double> const& coeffs)
{
    EXPECT_EQ(1u, coeffs.size());
    auto const c = coeffs[0];
    return [c](std::array<double, 3> const& /*coords*/) { return c; };
}

LocalAssemblerDataInterface::Function getFLin(std::vector<double> const& coeffs)
{
    EXPECT_EQ(4u, coeffs.size());
    return [coeffs](std::array<double, 3> const& coords) {
        auto const x = coords[0];
        auto const y = coords[1];
        auto const z = coords[2];
        return coeffs[0] + coeffs[1] * x + coeffs[2] * y + coeffs[3] * z;
    };
}

LocalAssemblerDataInterface::Function getFQuad(
    std::vector<double> const& coeffs)
{
    EXPECT_EQ(10u, coeffs.size());
    return [coeffs](std::array<double, 3> const& coords) {
        auto const x = coords[0];
        auto const y = coords[1];
        auto const z = coords[2];
        return coeffs[0] + coeffs[1] * x + coeffs[2] * y + coeffs[3] * z +
               coeffs[4] * x * y + coeffs[5] * y * z + coeffs[6] * z * x +
               coeffs[7] * x * x + coeffs[8] * y * y + coeffs[9] * z * z;
    };
}

std::pair<std::vector<double>, LocalAssemblerDataInterface::Function> getF(
    unsigned polynomial_order)
{
    std::vector<double> coeffs;

    switch (polynomial_order)
    {
        case 0:
            coeffs.resize(1);
            fillVectorRandomly(coeffs);
            return {coeffs, getFConst(coeffs)};
        case 1:
            coeffs.resize(4);
            fillVectorRandomly(coeffs);
            return {coeffs, getFLin(coeffs)};
        case 2:
            coeffs.resize(10);
            fillVectorRandomly(coeffs);
            return {coeffs, getFQuad(coeffs)};
    }

    OGS_FATAL("unsupported polynomial order: %d.", polynomial_order);
}

TEST(MathLib, IntegrationGaussLegendreHexConst)
{
    auto const eps = 2 * std::numeric_limits<double>::epsilon();
    std::unique_ptr<MeshLib::Mesh> mesh_hex(
        MeshLib::IO::VtuInterface::readVTUFile(BaseLib::BuildInfo::data_path +
                                               "/MathLib/unit_cube_hex.vtu"));

    for (unsigned integration_order : {1, 2, 3})
    {
        DBUG("\n==== integration order: %u.\n", integration_order);
        TestProcess pcs_hex(*mesh_hex, integration_order);

        const unsigned polynomial_order = 0;
        auto f = getF(polynomial_order);

        auto const integral_hex = pcs_hex.integrate(f.second);
        EXPECT_NEAR(f.first[0], integral_hex, eps);
    }
}

TEST(MathLib, IntegrationGaussLegendreTetConst)
{
    auto const eps = std::numeric_limits<double>::epsilon();
    std::unique_ptr<MeshLib::Mesh> mesh_tet(
        MeshLib::IO::VtuInterface::readVTUFile(BaseLib::BuildInfo::data_path +
                                               "/MathLib/unit_cube_tet.vtu"));

    for (unsigned integration_order : {1, 2, 3})
    {
        DBUG("\n==== integration order: %u.\n", integration_order);
        TestProcess pcs_tet(*mesh_tet, integration_order);

        const unsigned polynomial_order = 0;
        auto f = getF(polynomial_order);

        auto const integral_tet = pcs_tet.integrate(f.second);
        EXPECT_NEAR(f.first[0], integral_tet, eps);
    }
}

TEST(MathLib, IntegrationGaussLegendreTet)
{
    auto const eps = std::numeric_limits<double>::epsilon();
    std::unique_ptr<MeshLib::Mesh> mesh_tet(
        MeshLib::IO::VtuInterface::readVTUFile(BaseLib::BuildInfo::data_path +
                                               "/MathLib/unit_cube_tet.vtu"));
    std::unique_ptr<MeshLib::Mesh> mesh_hex(
        MeshLib::IO::VtuInterface::readVTUFile(BaseLib::BuildInfo::data_path +
                                               "/MathLib/unit_cube_hex.vtu"));

    for (unsigned integration_order : {1, 2, 3})
    {
        DBUG("\n==== integration order: %u.\n", integration_order);
        TestProcess pcs_tet(*mesh_tet, integration_order);
        TestProcess pcs_hex(*mesh_hex, integration_order);

        for (unsigned polynomial_order : {0, 1, 2})
        {
            if (polynomial_order > 2 * integration_order - 1)
                break;

            DBUG("  == polynomial order: %u.", polynomial_order);
            auto f = getF(polynomial_order);

            auto const integral_tet = pcs_tet.integrate(f.second);
            auto const integral_hex = pcs_hex.integrate(f.second);
            // DBUG("  integrals: %g, %g.", integral_tet, integral_hex);
            EXPECT_NEAR(integral_hex, integral_tet, eps);
        }
    }
}
