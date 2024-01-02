/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifdef OGS_USE_MFRONT

#include <gtest/gtest.h>

#include <numeric>

#include "MaterialLib/MPL/Utils/Tensor.h"
#include "MaterialLib/SolidModels/MFront/TangentOperatorBlocksView.h"
#include "MaterialLib/SolidModels/MFront/Variable.h"
#include "OGSMFrontTestVariables.h"
#include "Tests/TestTools.h"

template <class Dim>
struct MaterialLib_TangentOperatorBlocksView : ::testing::Test
{
};

using MaterialLib_TangentOperatorBlocksView_TestCases =
    ::testing::Types<std::integral_constant<int, 2>,
                     std::integral_constant<int, 3>>;

TYPED_TEST_SUITE(MaterialLib_TangentOperatorBlocksView,
                 MaterialLib_TangentOperatorBlocksView_TestCases);

TYPED_TEST(MaterialLib_TangentOperatorBlocksView, Test1)
{
    namespace MB = mgis::behaviour;
    using Var = MB::Variable;
    using namespace boost::mp11;
    namespace MSM = MaterialLib::Solids::MFront;

    constexpr int dim = TypeParam::value;
    constexpr int kv_size =
        MathLib::KelvinVector::kelvin_vector_dimensions(dim);
    constexpr int tensor_size = MaterialPropertyLib::tensorSize(dim);

    const std::vector to_blocks{
        // dsigma/dtensor
        std::pair{Var{"Stress", Var::STENSOR}, Var{"tensor", Var::TENSOR}},
        // dvector/dp
        std::pair{Var{"vector", Var::VECTOR},
                  Var{"LiquidPressure", Var::SCALAR}},
        // dsigma/dT
        std::pair{Var{"Stress", Var::STENSOR},
                  Var{"Temperature", Var::SCALAR}}};

    const std::size_t total_data_size =
        kv_size * tensor_size + dim * 1 + kv_size * 1;

    using Gradients = mp_list<MSM::LiquidPressure, Tensor>;
    using TDynForces = mp_list<Vector, MSM::Stress>;
    using ExtStateVars = mp_list<MSM::Temperature>;

    MSM::OGSMFrontTangentOperatorBlocksView<
        dim,
        MSM::ForcesGradsCombinations<Gradients, TDynForces, ExtStateVars>::type>
        view(to_blocks);

    MSM::OGSMFrontTangentOperatorData data{};
    data.data.resize(total_data_size);
    std::iota(begin(data.data), end(data.data), 0);

    // dvector/dp != 0
    {
        using Vec = Eigen::Vector<double, dim>;
        auto const offset = kv_size * tensor_size;
        Vec const expected = Vec::LinSpaced(offset, offset + dim - 1);

        auto const b = view.block(Vector{}, MSM::liquid_pressure, data);

        EXPECT_PRED_FORMAT2(Tests::EigenIsNear{}, expected, b);
    }

    // dvector/dtensor = 0
    {
        using Mat = Eigen::Matrix<double, dim, tensor_size>;
        Mat const zero = Mat::Zero();

        auto const b = view.block(Vector{}, Tensor{}, data);

        EXPECT_PRED_FORMAT2(Tests::EigenIsNear{}, zero, b);
    }

    // dvector/T = 0
    {
        using Vec = Eigen::Vector<double, dim>;
        Vec const zero = Vec::Zero();

        auto const b = view.block(Vector{}, MSM::temperature, data);

        EXPECT_PRED_FORMAT2(Tests::EigenIsNear{}, zero, b);
    }

    // dsigma/dp = 0
    {
        using Vec = Eigen::Vector<double, kv_size>;
        Vec const zero = Vec::Zero();

        auto const b = view.block(MSM::stress, MSM::liquid_pressure, data);

        EXPECT_PRED_FORMAT2(Tests::EigenIsNear{}, zero, b);
    }

    // dsigma/dtensor != 0
    {
        using Mat = Eigen::Matrix<double, kv_size, tensor_size>;
        Mat expected = Mat::Zero();
        for (int r = 0; r < kv_size; ++r)
        {
            for (int c = 0; c < tensor_size; ++c)
            {
                expected(r, c) = tensor_size * r + c;
            }
        }

        auto const b = view.block(MSM::stress, Tensor{}, data);

        EXPECT_PRED_FORMAT2(Tests::EigenIsNear{}, expected, b);
    }

    // dsigma/dT != 0
    {
        using Vec = Eigen::Vector<double, kv_size>;
        auto const offset = kv_size * tensor_size + dim;
        Vec const expected = Vec::LinSpaced(offset, offset + kv_size - 1);

        auto const b = view.block(MSM::stress, MSM::temperature, data);

        EXPECT_PRED_FORMAT2(Tests::EigenIsNear{}, expected, b);
    }
}

TYPED_TEST(MaterialLib_TangentOperatorBlocksView, Test2)
{
    namespace MB = mgis::behaviour;
    using Var = MB::Variable;
    using namespace boost::mp11;
    namespace MSM = MaterialLib::Solids::MFront;

    constexpr int dim = TypeParam::value;
    constexpr int kv_size =
        MathLib::KelvinVector::kelvin_vector_dimensions(dim);
    constexpr int tensor_size = MaterialPropertyLib::tensorSize(dim);

    const std::vector to_blocks{
        // dsigma/dtensor
        std::pair{Var{"Stress", Var::STENSOR}, Var{"tensor", Var::TENSOR}},
        // dvector/dp
        std::pair{Var{"Saturation", Var::VECTOR},
                  Var{"LiquidPressure", Var::SCALAR}}};

    const std::size_t total_data_size = kv_size * tensor_size + dim * 1;

    using Gradients = mp_list<MSM::LiquidPressure, Tensor>;
    using TDynForces = mp_list<MSM::Saturation, MSM::Stress>;
    using ExtStateVars = mp_list<MSM::Temperature>;

    MSM::OGSMFrontTangentOperatorBlocksView<
        dim,
        MSM::ForcesGradsCombinations<Gradients, TDynForces, ExtStateVars>::type>
        view(to_blocks);

    MSM::OGSMFrontTangentOperatorData data{};
    data.data.resize(total_data_size);
    std::iota(begin(data.data), end(data.data), 0);

    // dsat/dp != 0
    {
        auto const offset = kv_size * tensor_size;
        auto const expected = offset;

        auto const b = view.block(MSM::saturation, MSM::liquid_pressure, data);

        // double comparison (not Eigen::Map), not contained in the first test
        // case
        EXPECT_DOUBLE_EQ(expected, b);
    }

    // dsat/dtensor = 0
    {
        using Vec = Eigen::RowVector<double, tensor_size>;
        Vec const zero = Vec::Zero();

        auto const b = view.block(MSM::saturation, Tensor{}, data);

        EXPECT_PRED_FORMAT2(Tests::EigenIsNear{}, zero, b);
    }

    // dsat/dT = 0
    {
        auto const b = view.block(MSM::saturation, MSM::temperature, data);

        EXPECT_DOUBLE_EQ(0, b);
    }

    // dsigma/dp = 0
    {
        using Vec = Eigen::Vector<double, kv_size>;
        Vec const zero = Vec::Zero();

        auto const b = view.block(MSM::stress, MSM::liquid_pressure, data);

        EXPECT_PRED_FORMAT2(Tests::EigenIsNear{}, zero, b);
    }

    // dsigma/dtensor != 0
    {
        using Mat = Eigen::Matrix<double, kv_size, tensor_size>;
        Mat expected = Mat::Zero();
        for (int r = 0; r < kv_size; ++r)
        {
            for (int c = 0; c < tensor_size; ++c)
            {
                expected(r, c) = tensor_size * r + c;
            }
        }

        auto const b = view.block(MSM::stress, Tensor{}, data);

        EXPECT_PRED_FORMAT2(Tests::EigenIsNear{}, expected, b);
    }

    // dsigma/dT == 0
    {
        using Vec = Eigen::Vector<double, kv_size>;
        Vec const zero = Vec::Zero();

        auto const b = view.block(MSM::stress, MSM::temperature, data);

        EXPECT_PRED_FORMAT2(Tests::EigenIsNear{}, zero, b);
    }
}

TYPED_TEST(MaterialLib_TangentOperatorBlocksView, UnusedMFrontBlock)
{
    namespace MB = mgis::behaviour;
    using Var = MB::Variable;
    using namespace boost::mp11;
    namespace MSM = MaterialLib::Solids::MFront;

    constexpr int dim = TypeParam::value;

    const std::vector to_blocks{
        // dsigma/dtensor
        std::pair{Var{"Stress", Var::STENSOR}, Var{"tensor", Var::TENSOR}},
        // dvector/dp
        std::pair{Var{"Saturation", Var::VECTOR},
                  Var{"LiquidPressure", Var::SCALAR}}};

    using Gradients = mp_list<MSM::LiquidPressure>;
    using TDynForces = mp_list<MSM::Saturation, MSM::Stress>;
    using ExtStateVars = mp_list<MSM::Temperature>;

    using View = MSM::OGSMFrontTangentOperatorBlocksView<
        dim,
        MSM::ForcesGradsCombinations<Gradients, TDynForces, ExtStateVars>::
            type>;
    ASSERT_ANY_THROW({ View view(to_blocks); });
    // The constructor of View calls OGS_FATAL if MFront provides a tangent
    // operator block for which there is no entry in ForcesGradsCombinations,
    // i.e., that cannot be accessed in OGS through the view.
    // In this test dStress/dtensor is such an inaccessible block.
}

TYPED_TEST(MaterialLib_TangentOperatorBlocksView, OverspecifiedMFrontBlock)
{
    namespace MB = mgis::behaviour;
    using Var = MB::Variable;
    using namespace boost::mp11;
    namespace MSM = MaterialLib::Solids::MFront;

    constexpr int dim = TypeParam::value;

    const std::vector to_blocks{
        // dsigma/dtensor
        std::pair{Var{"Stress", Var::STENSOR}, Var{"tensor", Var::TENSOR}},
        std::pair{Var{"Stress", Var::STENSOR}, Var{"Temperature", Var::SCALAR}},
        // dvector/dp
        std::pair{Var{"Saturation", Var::VECTOR},
                  Var{"LiquidPressure", Var::SCALAR}}};

    using Gradients = mp_list<MSM::LiquidPressure, Tensor, Vector>;
    using TDynForces = mp_list<MSM::Saturation, MSM::Stress>;
    using ExtStateVars = mp_list<MSM::Temperature>;

    using View = MSM::OGSMFrontTangentOperatorBlocksView<
        dim,
        MSM::ForcesGradsCombinations<Gradients, TDynForces, ExtStateVars>::
            type>;
    ASSERT_NO_THROW(View{to_blocks});
    // to_blocks covers only a subset of all of the ForcesGradsCombinations, but
    // that's OK, since not all MFront models might provide all blocks.
    ASSERT_NO_THROW(View{to_blocks});
}
#endif
