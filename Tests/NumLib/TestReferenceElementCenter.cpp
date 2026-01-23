
#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <boost/mp11.hpp>

#include "NumLib/Fem/FiniteElement/ElementTraitsLagrange.h"
#include "NumLib/Fem/InitShapeMatrices.h"
#include "NumLib/Fem/Integration/IntegrationMethodRegistry.h"
#include "NumLib/Fem/ReferenceElement.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"

template <typename ShapeFunction>
struct ReferenceElementCenterTest : public ::testing::Test
{
};

template <typename ElementTraitsLagrange>
using GetShapeFunctionFromElementTraitsLagrange =
    typename ElementTraitsLagrange::ShapeFunction;

using ReferenceElementCenterTestTypes = boost::mp11::mp_rename<
    boost::mp11::mp_transform<GetShapeFunctionFromElementTraitsLagrange,
                              NumLib::AllElementTraitsLagrange>,
    testing::Types>;

TYPED_TEST_SUITE(ReferenceElementCenterTest, ReferenceElementCenterTestTypes);

template <typename ShapeFunction, typename MeshElementType>
Eigen::Vector3d centerOfGravityByIntegral(MeshElementType const& element)
{
    constexpr int dim = ShapeFunction::DIM > 0 ? ShapeFunction::DIM : 1;

    using ShapeMatricesType = ShapeMatrixPolicyType<ShapeFunction, dim>;

    auto const& int_meth =
        NumLib::IntegrationMethodRegistry::getIntegrationMethod<
            MeshElementType>(NumLib::IntegrationOrder{4});

    auto const Ns = NumLib::initShapeMatrices<ShapeFunction,
                                              ShapeMatricesType,
                                              dim,
                                              NumLib::ShapeMatrixType::N_J>(
        element, false /*is_axially_symmetric*/, int_meth);

    auto r_G = Eigen::Vector3d::Zero().eval();
    double V = 0;

    unsigned const n_ips = int_meth.getNumberOfPoints();
    for (unsigned ip = 0; ip < n_ips; ++ip)
    {
        auto const& wp = int_meth.getWeightedPoint(ip);
        auto const& N_J = Ns[ip];
        auto const& N = N_J.N;
        double const detJ = N_J.detJ;
        double const w = wp.getWeight();

        auto const r =
            NumLib::interpolateCoordinates<ShapeFunction, ShapeMatricesType>(
                element, N);

        V += detJ * w;
        r_G += Eigen::Map<Eigen::Vector3d const>(r.data()) * detJ * w;
    }

    r_G /= V;

    return r_G;
}

TYPED_TEST(ReferenceElementCenterTest, ReferenceElementCenter)
{
    using ShapeFunction = TypeParam;
    using MeshElementType = ShapeFunction::MeshElement;

    constexpr int dim = ShapeFunction::DIM > 0 ? ShapeFunction::DIM : 1;

    using ShapeMatricesType = ShapeMatrixPolicyType<ShapeFunction, dim>;

    NumLib::ReferenceElement<MeshElementType> reference_element;
    auto const& element = reference_element.element;

    auto const N_centre =
        NumLib::initShapeMatricesAtElementCenter<ShapeFunction,
                                                 ShapeMatricesType,
                                                 dim,
                                                 NumLib::ShapeMatrixType::N>(
            element, false /* is_axially_symmetric */)
            .N;

    auto const center_actual =
        NumLib::interpolateCoordinates<ShapeFunction, ShapeMatricesType>(
            element, N_centre);

    auto const center_expected =
        centerOfGravityByIntegral<ShapeFunction>(element);

    EXPECT_THAT(center_actual,
                testing::Pointwise(
                    testing::DoubleNear(std::numeric_limits<double>::epsilon()),
                    center_expected));
}
