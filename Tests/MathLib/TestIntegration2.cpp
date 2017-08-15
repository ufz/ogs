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
            // We need to assert that detJ is constant in each element in order
            // to be sure that in every element that we are integrating over the
            // effective polynomial degree is the one we expect.
            EXPECT_NEAR(sms[0].detJ, sms[ip].detJ,
                        std::numeric_limits<double>::epsilon());
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

struct FBase
{
private:
    static std::vector<double> initCoeffs(std::size_t const num_coeffs)
    {
        std::vector<double> cs(num_coeffs);
        fillVectorRandomly(cs);
        return cs;
    }

public:
    FBase(std::size_t const num_coeffs) : coeffs(initCoeffs(num_coeffs)) {}

    virtual double operator()(
        std::array<double, 3> const& /*coords*/) const = 0;
    virtual double getAnalyticalIntegralOverUnitCube() const = 0;

    LocalAssemblerDataInterface::Function getClosure() const
    {
        return [this](std::array<double, 3> const& coords) {
            return this->operator()(coords);
        };
    }

    virtual ~FBase() = default;

    std::vector<double> const coeffs;
};

struct FConst final : FBase
{
    FConst() : FBase(1) {}

    double operator()(std::array<double, 3> const&) const override
    {
        return coeffs[0];
    }

    double getAnalyticalIntegralOverUnitCube() const override
    {
        return coeffs[0];
    }
};

struct FLin final : FBase
{
    FLin() : FBase(4) {}

    double operator()(std::array<double, 3> const& coords) const override
    {
        auto const x = coords[0];
        auto const y = coords[1];
        auto const z = coords[2];
        return coeffs[0] + coeffs[1] * x + coeffs[2] * y + coeffs[3] * z;
    }

    double getAnalyticalIntegralOverUnitCube() const override
    {
        return coeffs[0];
    }
};

struct FQuad final : FBase
{
    FQuad() : FBase(10) {}

    double operator()(std::array<double, 3> const& coords) const override
    {
        auto const x = coords[0];
        auto const y = coords[1];
        auto const z = coords[2];
        return coeffs[0] + coeffs[1] * x + coeffs[2] * y + coeffs[3] * z +
               coeffs[4] * x * y + coeffs[5] * y * z + coeffs[6] * z * x +
               coeffs[7] * x * x + coeffs[8] * y * y + coeffs[9] * z * z;
    }

    double getAnalyticalIntegralOverUnitCube() const override
    {
        double const a = -.5;
        double const b = .5;

        double const a3 = a * a * a;
        double const b3 = b * b * b;

        return coeffs[0] + (coeffs[7] + coeffs[8] + coeffs[9]) * (b3 - a3) / 3.;
    }
};

struct F3DSeparablePolynomial final : FBase
{
    F3DSeparablePolynomial(unsigned polynomial_degree)
        : FBase(3 * polynomial_degree + 3), _degree(polynomial_degree)
    {
    }

    // f(x, y, z) = g(x) * h(y) * i(z)
    double operator()(std::array<double, 3> const& coords) const override
    {
        double res = 1.0;
        for (unsigned d : {0, 1, 2})
        {
            auto const x = coords[d];

            double poly = 0.0;
            for (unsigned n = 0; n <= _degree; ++n)
            {
                poly += coeffs[n + d * (_degree + 1)] * std::pow(x, n);
            }

            res *= poly;
        }

        return res;
    }

    // [ F(x, y, z) ]_a^b = [ G(x) ]_a^b * [ H(y) ]_a^b * [ I(z) ]_a^b
    double getAnalyticalIntegralOverUnitCube() const override
    {
        double const a = -.5;
        double const b = .5;

        double res = 1.0;
        for (unsigned d : {0, 1, 2})
        {
            double poly = 0.0;
            for (unsigned n = 0; n <= _degree; ++n)
            {
                poly += coeffs[n + d * (_degree + 1)] *
                        (std::pow(b, n + 1) - std::pow(a, n + 1)) / (n + 1);
            }

            res *= poly;
        }

        return res;
    }

private:
    unsigned const _degree;
};

std::unique_ptr<FBase> getF(unsigned polynomial_order)
{
    std::vector<double> coeffs;

    switch (polynomial_order)
    {
        case 0:
            return std::make_unique<FConst>();
        case 1:
            return std::make_unique<FLin>();
        case 2:
            return std::make_unique<FQuad>();
    }

    OGS_FATAL("unsupported polynomial order: %d.", polynomial_order);
}

TEST(MathLib, IntegrationGaussLegendreTet)
{
    auto const eps = 2 * std::numeric_limits<double>::epsilon();
    std::unique_ptr<MeshLib::Mesh> mesh_tet(
        MeshLib::IO::VtuInterface::readVTUFile(BaseLib::BuildInfo::data_path +
                                               "/MathLib/unit_cube_tet.vtu"));

    for (unsigned integration_order : {1, 2, 3})
    {
        DBUG("\n==== integration order: %u.\n", integration_order);
        TestProcess pcs_tet(*mesh_tet, integration_order);

        for (unsigned polynomial_order : {0, 1, 2})
        {
            if (polynomial_order > 2 * integration_order - 1)
                break;

            DBUG("  == polynomial order: %u.", polynomial_order);
            auto f = getF(polynomial_order);

            auto const integral_tet = pcs_tet.integrate(f->getClosure());
            EXPECT_NEAR(f->getAnalyticalIntegralOverUnitCube(), integral_tet,
                        eps);
        }
    }
}

TEST(MathLib, IntegrationGaussLegendreHex)
{
    auto const eps = 2 * std::numeric_limits<double>::epsilon();
    std::unique_ptr<MeshLib::Mesh> mesh_hex(
        MeshLib::IO::VtuInterface::readVTUFile(BaseLib::BuildInfo::data_path +
                                               "/MathLib/unit_cube_hex.vtu"));

    for (unsigned integration_order : {1, 2, 3})
    {
        DBUG("\n==== integration order: %u.\n", integration_order);
        TestProcess pcs_hex(*mesh_hex, integration_order);

        for (unsigned polynomial_order : {0, 1, 2})
        {
            if (polynomial_order > 2 * integration_order - 1)
                break;

            DBUG("  == polynomial order: %u.", polynomial_order);
            auto f = getF(polynomial_order);

            auto const integral_hex = pcs_hex.integrate(f->getClosure());
            EXPECT_NEAR(f->getAnalyticalIntegralOverUnitCube(), integral_hex,
                        eps);
        }
    }
}

TEST(MathLib, IntegrationGaussLegendreTetSeparablePolynomial)
{
    auto const eps = 2 * std::numeric_limits<double>::epsilon();

    std::unique_ptr<MeshLib::Mesh> mesh_tet(
        MeshLib::IO::VtuInterface::readVTUFile(BaseLib::BuildInfo::data_path +
                                               "/MathLib/unit_cube_tet.vtu"));

    for (unsigned integration_order : {1, 2, 3})
    {
        DBUG("\n==== integration order: %u.\n", integration_order);
        TestProcess pcs_tet(*mesh_tet, integration_order);

        for (unsigned polynomial_order = 0;
             // Gauss-Legendre integration is exact up to this order!
             polynomial_order < 2 * integration_order;
             ++polynomial_order)
        {
            DBUG("  == polynomial order: %u.", polynomial_order);
            F3DSeparablePolynomial f(polynomial_order);

            auto const integral_tet = pcs_tet.integrate(f.getClosure());
            EXPECT_NEAR(f.getAnalyticalIntegralOverUnitCube(), integral_tet,
                        eps);
        }
    }
}

TEST(MathLib, IntegrationGaussLegendreHexSeparablePolynomial)
{
    auto const eps = 2 * std::numeric_limits<double>::epsilon();

    std::unique_ptr<MeshLib::Mesh> mesh_hex(
        MeshLib::IO::VtuInterface::readVTUFile(BaseLib::BuildInfo::data_path +
                                               "/MathLib/unit_cube_hex.vtu"));

    for (unsigned integration_order : {1, 2, 3, 4})
    {
        DBUG("\n==== integration order: %u.\n", integration_order);
        TestProcess pcs_hex(*mesh_hex, integration_order);

        for (unsigned polynomial_order = 0;
             // Gauss-Legendre integration is exact up to this order!
             polynomial_order < 2 * integration_order;
             ++polynomial_order)
        {
            DBUG("  == polynomial order: %u.", polynomial_order);
            F3DSeparablePolynomial f(polynomial_order);

            auto const integral_hex = pcs_hex.integrate(f.getClosure());
            EXPECT_NEAR(f.getAnalyticalIntegralOverUnitCube(), integral_hex,
                        eps);
        }
    }
}
