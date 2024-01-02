/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "MaterialLib/MPL/Utils/Tensor.h"
#include "MaterialLib/SolidModels/MFront/ThermodynamicForcesView.h"
#include "MathLib/KelvinVector.h"
#include "OGSMFrontTestVariables.h"
#include "Tests/TestTools.h"

namespace MSM = MaterialLib::Solids::MFront;

// /////////////////////////////////////////////////////////////////////////////
// The unit tests in this file all have a version that needs MFront and a
// version that doesn't. For the latter we have to provide the right mock
// variables, here.
// The mocking does not only allow tests in environments that do not use MFront,
// but also shows that OGSMFrontThermodynamicForcesView works entirely without
// MFront, i.e., can in principle be used in other contexts, too. The necessary
// adapters have to be like the mock variables in this file.

struct MockStress
{
    template <int DisplacementDim>
    static constexpr std::size_t size()
    {
        return rows<DisplacementDim>() * cols<DisplacementDim>();
    }

    template <int DisplacementDim>
    static constexpr std::size_t rows()
    {
        return MathLib::KelvinVector::kelvin_vector_dimensions(DisplacementDim);
    }

    template <int DisplacementDim>
    static constexpr std::size_t cols()
    {
        return 1;
    }
};

struct MockSaturation
{
    template <int DisplacementDim>
    static constexpr std::size_t size()
    {
        return rows<DisplacementDim>() * cols<DisplacementDim>();
    }

    template <int DisplacementDim>
    static constexpr std::size_t rows()
    {
        return 1;
    }

    template <int DisplacementDim>
    static constexpr std::size_t cols()
    {
        return 1;
    }
};

struct MockVector
{
    template <int DisplacementDim>
    static constexpr std::size_t size()
    {
        return rows<DisplacementDim>() * cols<DisplacementDim>();
    }

    template <int DisplacementDim>
    static constexpr std::size_t rows()
    {
        return DisplacementDim;
    }

    template <int DisplacementDim>
    static constexpr std::size_t cols()
    {
        return 1;
    }
};

struct MockTensor
{
    template <int DisplacementDim>
    static constexpr std::size_t size()
    {
        return rows<DisplacementDim>() * cols<DisplacementDim>();
    }

    template <int DisplacementDim>
    static constexpr std::size_t rows()
    {
        return MaterialPropertyLib::tensorSize(DisplacementDim);
    }

    template <int DisplacementDim>
    static constexpr std::size_t cols()
    {
        return 1;
    }
};

template <typename Stress, typename Saturation>
static void
test_MaterialLib_ThermodynamicForcesView_ReadAccess_STensorScalar_3D_impl()
{
    using TDynForces = boost::mp11::mp_list<Stress, Saturation>;

    MSM::OGSMFrontThermodynamicForcesData const data{{
        1, 2, 3, 4, 5, 6,  // stress
        7                  // saturation
    }};

    MSM::OGSMFrontThermodynamicForcesView<3, TDynForces> view;

    auto const stress = view.block(Stress{}, data);

    EXPECT_THAT(stress,
                testing::Pointwise(testing::DoubleEq(), {1, 2, 3, 4, 5, 6}));

    auto const saturation = view.block(Saturation{}, data);

    EXPECT_DOUBLE_EQ(7, saturation);
}

#ifdef OGS_USE_MFRONT
TEST(MaterialLib_ThermodynamicForcesView, ReadAccess_STensorScalar_3D)
{
    test_MaterialLib_ThermodynamicForcesView_ReadAccess_STensorScalar_3D_impl<
        MSM::Stress,
        MSM::Saturation>();
}
#endif

TEST(MaterialLib_ThermodynamicForcesView, ReadAccess_STensorScalar_3D_Mock)
{
    test_MaterialLib_ThermodynamicForcesView_ReadAccess_STensorScalar_3D_impl<
        MockStress,
        MockSaturation>();
}

template <typename Vector, typename Tensor>
static void
test_MaterialLib_ThermodynamicForcesView_ReadAccess_VectorTensor_2D_impl()
{
    using TDynForces = boost::mp11::mp_list<Vector, Tensor>;

    MSM::OGSMFrontThermodynamicForcesData const data{{
        1, 2,          // vector data
        3, 4, 5, 6, 7  // tensor data
    }};

    MSM::OGSMFrontThermodynamicForcesView<2, TDynForces> view;

    auto const vector = view.block(Vector{}, data);

    EXPECT_THAT(vector, testing::Pointwise(testing::DoubleEq(), {1, 2}));

    auto const tensor = view.block(Tensor{}, data);
    auto const tensor_expected =
        (Eigen::Matrix<double, 5, 1>{} << 3, 4, 5, 6, 7).finished();

    ASSERT_PRED_FORMAT2(Tests::EigenIsNear{}, tensor, tensor_expected);
}

#ifdef OGS_USE_MFRONT
TEST(MaterialLib_ThermodynamicForcesView, ReadAccess_VectorTensor_2D)
{
    test_MaterialLib_ThermodynamicForcesView_ReadAccess_VectorTensor_2D_impl<
        Vector,
        Tensor>();
}
#endif

TEST(MaterialLib_ThermodynamicForcesView, ReadAccess_VectorTensor_2D_Mock)
{
    test_MaterialLib_ThermodynamicForcesView_ReadAccess_VectorTensor_2D_impl<
        MockVector,
        MockTensor>();
}

template <typename Saturation, typename Vector>
static void
test_MaterialLib_ThermodynamicForcesView_WriteAccess_ScalarVector_2D_impl()
{
    using TDynForces = boost::mp11::mp_list<Saturation, Vector>;

    MSM::OGSMFrontThermodynamicForcesData data{{
        1,    // saturation
        2, 3  // vector data
    }};

    MSM::OGSMFrontThermodynamicForcesView<2, TDynForces> view;

    view.block(Saturation{}, data) = 5;
    view.block(Vector{}, data)[1] = 7;

    EXPECT_THAT(data.data, testing::Pointwise(testing::DoubleEq(), {5, 2, 7}));
}

#ifdef OGS_USE_MFRONT
TEST(MaterialLib_ThermodynamicForcesView, WriteAccess_ScalarVector_2D)
{
    test_MaterialLib_ThermodynamicForcesView_WriteAccess_ScalarVector_2D_impl<
        MSM::Saturation,
        Vector>();
}
#endif

TEST(MaterialLib_ThermodynamicForcesView, WriteAccess_ScalarVector_2D_Mock)
{
    test_MaterialLib_ThermodynamicForcesView_WriteAccess_ScalarVector_2D_impl<
        MockSaturation,
        MockVector>();
}

template <typename Tensor, typename STensor>
static void
test_MaterialLib_ThermodynamicForcesView_WriteAccess_TensorSTensor_3D_impl()
{
    using TDynForces = boost::mp11::mp_list<Tensor, STensor>;

    MSM::OGSMFrontThermodynamicForcesData data{{
        1, 2, 3, 4, 5, 6, 7, 8, 9,  // tensor data
        10, 11, 12, 13, 14, 15      // strain
    }};

    MSM::OGSMFrontThermodynamicForcesView<3, TDynForces> view;

    view.block(STensor{}, data).template segment<3>(1) =
        Eigen::Vector3d(21, 22, 23);
    view.block(Tensor{}, data).template segment<2>(1) =
        Eigen::RowVector2d(31, 32);
    view.block(Tensor{}, data).template segment<2>(4) =
        Eigen::RowVector2d(41, 42);

    EXPECT_THAT(data.data,
                testing::Pointwise(
                    testing::DoubleEq(),
                    {1, 31, 32, 4, 41, 42, 7, 8, 9, 10, 21, 22, 23, 14, 15}));
}

#ifdef OGS_USE_MFRONT
TEST(MaterialLib_ThermodynamicForcesView, WriteAccess_TensorSTensor_3D)
{
    test_MaterialLib_ThermodynamicForcesView_WriteAccess_TensorSTensor_3D_impl<
        Tensor,
        MSM::Strain>();
}
#endif

TEST(MaterialLib_ThermodynamicForcesView, WriteAccess_TensorSTensor_3D_Mock)
{
    test_MaterialLib_ThermodynamicForcesView_WriteAccess_TensorSTensor_3D_impl<
        MockTensor,
        MockStress>();
}

template <typename Scalar, typename Vector, typename STensor, typename Tensor>
static void test_MaterialLib_ThermodynamicForcesView_DataSizes_impl()
{
    {
        using TDynForces = boost::mp11::mp_list<Tensor, STensor>;

        static_assert(MSM::OGSMFrontThermodynamicForcesView<3, TDynForces>::
                          data_size_all_forces == 9 + 6);

        static_assert(MSM::OGSMFrontThermodynamicForcesView<2, TDynForces>::
                          data_size_all_forces == 5 + 4);
    }

    // same as above, order changed
    {
        using TDynForces = boost::mp11::mp_list<STensor, Tensor>;

        static_assert(MSM::OGSMFrontThermodynamicForcesView<3, TDynForces>::
                          data_size_all_forces == 6 + 9);

        static_assert(MSM::OGSMFrontThermodynamicForcesView<2, TDynForces>::
                          data_size_all_forces == 4 + 5);
    }

    {
        using TDynForces = boost::mp11::mp_list<Scalar, Vector>;

        static_assert(MSM::OGSMFrontThermodynamicForcesView<3, TDynForces>::
                          data_size_all_forces == 1 + 3);

        static_assert(MSM::OGSMFrontThermodynamicForcesView<2, TDynForces>::
                          data_size_all_forces == 1 + 2);
    }

    // same as above, order changed
    {
        using TDynForces = boost::mp11::mp_list<Vector, Scalar>;

        static_assert(MSM::OGSMFrontThermodynamicForcesView<3, TDynForces>::
                          data_size_all_forces == 1 + 3);

        static_assert(MSM::OGSMFrontThermodynamicForcesView<2, TDynForces>::
                          data_size_all_forces == 1 + 2);
    }

    {
        using TDynForces = boost::mp11::mp_list<STensor, Scalar>;

        static_assert(MSM::OGSMFrontThermodynamicForcesView<3, TDynForces>::
                          data_size_all_forces == 6 + 1);

        static_assert(MSM::OGSMFrontThermodynamicForcesView<2, TDynForces>::
                          data_size_all_forces == 4 + 1);
    }

    // same as above, order changed
    {
        using TDynForces = boost::mp11::mp_list<Scalar, STensor>;

        static_assert(MSM::OGSMFrontThermodynamicForcesView<3, TDynForces>::
                          data_size_all_forces == 6 + 1);

        static_assert(MSM::OGSMFrontThermodynamicForcesView<2, TDynForces>::
                          data_size_all_forces == 4 + 1);
    }
}

#ifdef OGS_USE_MFRONT
TEST(MaterialLib_ThermodynamicForcesView, DataSizes)
{
    test_MaterialLib_ThermodynamicForcesView_DataSizes_impl<MSM::Saturation,
                                                            Vector,
                                                            MSM::Strain,
                                                            Tensor>();
}
#endif

TEST(MaterialLib_ThermodynamicForcesView, DataSizes_Mock)
{
    test_MaterialLib_ThermodynamicForcesView_DataSizes_impl<MockSaturation,
                                                            MockVector,
                                                            MockStress,
                                                            MockTensor>();
}
