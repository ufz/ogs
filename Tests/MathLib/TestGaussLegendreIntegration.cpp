/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <limits>

#include "InfoLib/TestInfo.h"
#include "MathLib/Integration/GaussLegendreTet.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/IO/VtkIO/VtuInterface.h"
#include "MeshLib/MeshSubset.h"

#ifdef USE_PETSC
#include "MeshLib/NodePartitionedMesh.h"
#endif

#include "NumLib/DOF/LocalToGlobalIndexMap.h"
#include "NumLib/Fem/InitShapeMatrices.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "Polynomials.h"
#include "ProcessLib/Utils/CreateLocalAssemblers.h"
#include "Tests/MeshLib/UnitCube.h"
#include "Tests/VectorUtils.h"

namespace GaussLegendreTest
{
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

class LocalAssemblerInterface
{
public:
    using Function = std::function<double(std::array<double, 3> const&)>;

    virtual double integrate(Function const& f,
                             unsigned const integration_order) const = 0;

    virtual ~LocalAssemblerInterface() = default;
};

template <typename ShapeFunction, typename IntegrationMethod, int GlobalDim>
class LocalAssembler : public LocalAssemblerInterface
{
    using ShapeMatricesType = ShapeMatrixPolicyType<ShapeFunction, GlobalDim>;
    using ShapeMatrices = typename ShapeMatricesType::ShapeMatrices;

public:
    LocalAssembler(MeshLib::Element const& e,
                   std::size_t const /*local_matrix_size*/,
                   bool is_axially_symmetric,
                   unsigned const /*integration_order*/)
        : _e(e)
    {
        if (is_axially_symmetric)
        {
            OGS_FATAL("Only testing Cartesian meshes!");
        }
    }

    double integrate(const Function& f,
                     unsigned const integration_order) const override
    {
        double integral = 0;

        auto const sms =
            NumLib::initShapeMatrices<ShapeFunction, ShapeMatricesType,
                                      GlobalDim>(
                _e, false /*is_axially_symmetric*/,
                IntegrationMethod{integration_order});
        IntegrationMethod integration_method{integration_order};

        for (unsigned ip = 0; ip < sms.size(); ++ip)
        {
            auto const& N = sms[ip].N;
            auto const weight =
                integration_method.getWeightedPoint(ip).getWeight();
            auto const coords = interpolateNodeCoordinates(_e, N);
            auto const function_value = f(coords);
            integral += function_value * sms[ip].detJ * weight;
        }

        return integral;
    }

private:
    MeshLib::Element const& _e;
};

class IntegrationTestProcess
{
public:
    IntegrationTestProcess(MeshLib::Mesh const& mesh,
                           unsigned const integration_order)
        : _integration_order(integration_order),
          _mesh_subset_all_nodes(mesh, mesh.getNodes())
    {
        std::vector<MeshLib::MeshSubset> all_mesh_subsets{
            _mesh_subset_all_nodes};

        _dof_table = std::make_unique<NumLib::LocalToGlobalIndexMap>(
            std::move(all_mesh_subsets), NumLib::ComponentOrder::BY_LOCATION);

        ProcessLib::createLocalAssemblers<LocalAssembler>(
            mesh.getDimension(), mesh.getElements(), *_dof_table,
            _local_assemblers, mesh.isAxiallySymmetric(), _integration_order);
    }

    double integrate(LocalAssemblerInterface::Function const& f) const
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

    std::vector<std::unique_ptr<LocalAssemblerInterface>> _local_assemblers;
};

}  // namespace GaussLegendreTest

/* *****************************************************************************
 *
 * The idea behind the tests in this file is to integrate polynomials of
 * different degree over the unit cube.
 *
 * Gauss-Legendre integration should be able to exactly integrate those up to a
 * certain degree.
 *
 * The coefficients of the tested polynomials are chosen randomly.
 *
 **************************************************************************** */

// Test fixture.
class MathLibIntegrationGaussLegendreTestBase : public ::testing::Test
{
protected:
    // Up to the polynomial order returned by this function the numerical
    // integration should be exact in the base case test.
    static unsigned maxExactPolynomialOrderBase(
        unsigned const integration_order)
    {
        // Base case is defined only for polynomials up to order 2.
        return std::min(2u, 2 * integration_order - 1);
    }

    // Up to the polynomial order returned by this function the numerical
    // integration should be exact in the separable polynomial test.
    static unsigned maxExactPolynomialOrderSeparable(
        unsigned const integration_order)
    {
        return 2 * integration_order - 1;
    }

    // Up to the polynomial order returned by this function the numerical
    // integration should be exact in the nonseparable polynomial test.
    static unsigned maxExactPolynomialOrderNonseparable(
        unsigned const integration_order)
    {
        return 2 * integration_order - 1;
    }
};

template <typename MeshElementType>
class MathLibIntegrationGaussLegendreTest
    : public MathLibIntegrationGaussLegendreTestBase
{
};

// Specialization for triangles.
template <>
class MathLibIntegrationGaussLegendreTest<MeshLib::Tri>
    : public MathLibIntegrationGaussLegendreTestBase
{
protected:
    static unsigned maxExactPolynomialOrderSeparable(
        unsigned const integration_order)
    {
        switch (integration_order)
        {
            case 1:
                return 0;  // constants
            case 2:
                return 1;  // most complicated monomial is x*y*z
            case 3:
                return 1;  // most complicated monomial is x*y*z
            case 4:
                return 2;  // most complicated monomial is x^2 * y^2 * z^2
        }

        OGS_FATAL("Unsupported integration order {}", integration_order);
    }

    static unsigned maxExactPolynomialOrderNonseparable(
        unsigned const integration_order)
    {
        switch (integration_order)
        {
            case 1:
                return 1;  // at most linear affine functions
            case 2:
                return 3;  // most complicated monomial is x^3 or x*y*z
            case 3:
                return 3;  // most complicated monomial is x^3 or x*y*z
            case 4:
                return 5;  // most complicated monomial is x^5 or x^2 * y^2 * z
        }

        OGS_FATAL("Unsupported integration order {}", integration_order);
    }
};

// Specialization for tetrahedra.
template <>
class MathLibIntegrationGaussLegendreTest<MeshLib::Tet>
    : public MathLibIntegrationGaussLegendreTestBase
{
protected:
    static unsigned maxExactPolynomialOrderSeparable(
        unsigned const integration_order)
    {
        switch (integration_order)
        {
            case 1:
                return 0;  // constants
            case 2:
                return 1;  // most complicated monomial is x*y*z
            case 3:
                return 1;  // most complicated monomial is x*y*z
            case 4:
                return 1;  // most complicated monomial is x*y*z
        }

        OGS_FATAL("Unsupported integration order {}", integration_order);
    }

    static unsigned maxExactPolynomialOrderNonseparable(
        unsigned const integration_order)
    {
        switch (integration_order)
        {
            case 1:
                return 1;  // at most linear affine functions
            case 2:
                return 3;  // most complicated monomial is x^3 or x*y*z
            case 3:
                return 5;  // most complicated monomial is x^5 or x^2 * y^2 * z
            case 4:
                return 5;  // most complicated monomial is x^5 or x^2 * y^2 * z
        }

        OGS_FATAL("Unsupported integration order {}", integration_order);
    }
};

// Specialization for prisms.
template <>
class MathLibIntegrationGaussLegendreTest<MeshLib::Prism>
    : public MathLibIntegrationGaussLegendreTestBase
{
protected:
    static unsigned maxExactPolynomialOrderSeparable(
        unsigned const integration_order)
    {
        switch (integration_order)
        {
            case 1:
                return 0;  // constants
            case 2:
                return 1;  // most complicated monomial is x*y*z
            case 3:
                return 2;  // most complicated monomial is x^2 * y^2 * z^2
            case 4:
                return 2;  // most complicated monomial is x^2 * y^2 * z^2
        }

        OGS_FATAL("Unsupported integration order {}", integration_order);
    }

    static unsigned maxExactPolynomialOrderNonseparable(
        unsigned const integration_order)
    {
        switch (integration_order)
        {
            case 1:
                return 1;  // at most linear affine functions
            case 2:
                return 3;  // most complicated monomial is x^3 or x*y*z
            case 3:
                return 5;  // most complicated monomial is x^5 or x^2 * y^2 * z
            case 4:
                return 5;  // most complicated monomial is x^5 or x^2 * y^2 * z
        }

        OGS_FATAL("Unsupported integration order {}", integration_order);
    }
};

// Specialization for prisms.
template <>
class MathLibIntegrationGaussLegendreTest<MeshLib::Pyramid>
    : public MathLibIntegrationGaussLegendreTestBase
{
protected:
    static unsigned maxExactPolynomialOrderSeparable(
        unsigned const integration_order)
    {
        // Note: In pyramids det J changes over the element. So effectively we
        // are integrating a polynomial of higher degree than indicated here.
        switch (integration_order)
        {
            case 1:
                return 1;  // most complicated monomial is x*y*z
            case 2:
                return 1;  // most complicated monomial is x*y*z
            case 3:
                return 3;  // most complicated monomial is x^3 * y^3 * z^3
            case 4:
                return 3;  // most complicated monomial is x^3 * y^3 * z^3
        }

        OGS_FATAL("Unsupported integration order {}", integration_order);
    }

    static unsigned maxExactPolynomialOrderNonseparable(
        unsigned const integration_order)
    {
        // Note: In pyramids det J changes over the element. So effectively we
        // are integrating a polynomial of higher degree than indicated here.
        switch (integration_order)
        {
            case 1:
                return 1;  // at most linear affine functions
            case 2:
                return 3;  // most complicated monomial is x^3 or x*y*z
            case 3:
                return 3;  // most complicated monomial is x^3 or x*y*z
            case 4:
                return 3;  // most complicated monomial is x^3 or x*y*z
        }

        OGS_FATAL("Unsupported integration order {}", integration_order);
    }
};

using MeshElementTypes =
    ::testing::Types<MeshLib::Line, MeshLib::Quad, MeshLib::Hex, MeshLib::Tri,
                     MeshLib::Tet, MeshLib::Prism, MeshLib::Pyramid>;

TYPED_TEST_SUITE(MathLibIntegrationGaussLegendreTest, MeshElementTypes);

template <typename MeshElementType, typename MaxExactPolynomialOrderFct,
          typename PolynomialFactory>
static void mathLibIntegrationGaussLegendreTestImpl(
    MaxExactPolynomialOrderFct const& max_exact_polynomial_order_fct,
    PolynomialFactory const& polynomial_factory)
{
    auto const eps = 10 * std::numeric_limits<double>::epsilon();
    auto const mesh_ptr = createUnitCube<MeshElementType>();

#ifndef USE_PETSC
    auto const& mesh = *mesh_ptr;
#else
    MeshLib::NodePartitionedMesh const mesh{*mesh_ptr};
#endif

    for (unsigned integration_order : {1, 2, 3, 4})
    {
        GaussLegendreTest::IntegrationTestProcess pcs(mesh, integration_order);

        for (unsigned polynomial_order = 0;
             polynomial_order <=
             max_exact_polynomial_order_fct(integration_order);
             ++polynomial_order)
        {
            auto const f = polynomial_factory(polynomial_order);

            auto const analytical_integral =
                f->getAnalyticalIntegralOverUnitCube();
            auto const numerical_integral = pcs.integrate(f->toFunction());
            EXPECT_NEAR(analytical_integral, numerical_integral, eps)
                << "integration order: " << integration_order
                << "\npolynomial order:  " << polynomial_order  //
                << "\nf:                 " << *f;
        }
    }
}

TYPED_TEST(MathLibIntegrationGaussLegendreTest, Base)
{
    using MeshElementType = TypeParam;
    auto const polynomial_factory = [](unsigned const polynomial_order)
    {
        auto constexpr dim = MeshElementType::dimension;
        return TestPolynomials::getF<dim>(polynomial_order);
    };

    mathLibIntegrationGaussLegendreTestImpl<MeshElementType>(
        this->maxExactPolynomialOrderBase, polynomial_factory);
}

TYPED_TEST(MathLibIntegrationGaussLegendreTest, SeparablePolynomial)
{
    using MeshElementType = TypeParam;

    auto const polynomial_factory = [](unsigned const polynomial_order)
    {
        auto constexpr dim = MeshElementType::dimension;
        return std::make_unique<TestPolynomials::FSeparablePolynomial<dim>>(
            polynomial_order);
    };

    mathLibIntegrationGaussLegendreTestImpl<MeshElementType>(
        this->maxExactPolynomialOrderSeparable, polynomial_factory);
}

TYPED_TEST(MathLibIntegrationGaussLegendreTest, NonseparablePolynomial)
{
    using MeshElementType = TypeParam;

    auto const polynomial_factory = [](unsigned const polynomial_order)
    {
        auto constexpr dim = MeshElementType::dimension;
        return std::make_unique<TestPolynomials::FNonseparablePolynomial<dim>>(
            polynomial_order);
    };

    mathLibIntegrationGaussLegendreTestImpl<MeshElementType>(
        this->maxExactPolynomialOrderNonseparable, polynomial_factory);
}
