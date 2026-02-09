// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "BaseLib/MPI.h"

#include <gtest/gtest.h>

using namespace BaseLib::MPI;

#ifdef USE_PETSC
// Fixture for MPI tests
struct MPI_BaseLib : public ::testing::Test
{
    Mpi mpi{};
};

TEST_F(MPI_BaseLib, CommunicatorRankAndSize)
{
    EXPECT_EQ(mpi.size, 3);
    EXPECT_GE(mpi.rank, 0);
    EXPECT_LT(mpi.rank, 3);
}

TEST_F(MPI_BaseLib, Allgather)
{
    std::vector<int> const gathered_values = allgather(mpi.rank, mpi);

    std::vector<int> const expected_values = {0, 1, 2};
    // Each rank's value should match its rank
    EXPECT_EQ(expected_values, gathered_values);
}

TEST_F(MPI_BaseLib, AllgatherVector)
{
    std::vector<int> const send_values(2, mpi.rank);
    std::vector<int> const gathered_values = allgather(send_values, mpi);

    std::vector<int> const expected_values = {0, 0, 1, 1, 2, 2};
    EXPECT_EQ(expected_values, gathered_values);
}

TEST_F(MPI_BaseLib, Allgatherv)
{
    int const rank_value =
        mpi.rank + 1;  // Each rank contributes rank+1 elements
    std::vector<int> const send_data(
        rank_value, mpi.rank);  // Each element is set to the rank number

    std::vector<int> gathered_values(1 + 2 + 3);

    std::vector<int> const offsets =
        allgatherv(std::span(send_data), gathered_values, mpi);
    std::vector<int> const expected_offsets = {0, 1, 3, 6};  // last element is
                                                             // size.
    std::vector<int> const expected_gathered_values = {0, 1, 1, 2, 2, 2};
    EXPECT_EQ(expected_offsets, offsets);
    EXPECT_EQ(expected_gathered_values, gathered_values);
}

TEST_F(MPI_BaseLib, Allreduce)
{
    int const sum = allreduce(mpi.rank + 1, MPI_SUM, mpi);
    int const expected_sum = 1 + 2 + 3;
    EXPECT_EQ(expected_sum, sum);

    int const product = allreduce(mpi.rank + 1, MPI_PROD, mpi);
    int const expected_product = 1 * 2 * 3;
    EXPECT_EQ(expected_product, product);

    int const min = allreduce(mpi.rank + 1, MPI_MIN, mpi);
    int const expected_min = 1;
    EXPECT_EQ(expected_min, min);
}

TEST_F(MPI_BaseLib, AllreduceVector)
{
    std::vector<int> const send_values = {mpi.rank + 1,
                                          mpi.size + mpi.rank + 1};
    std::vector<int> const sum = allreduce(send_values, MPI_SUM, mpi);
    std::vector<int> const expected_sum = {(1 + 2 + 3), (4 + 5 + 6)};
    EXPECT_EQ(expected_sum, sum);

    std::vector<int> const product = allreduce(send_values, MPI_PROD, mpi);
    std::vector<int> const expected_product = {(1 * 2 * 3), (4 * 5 * 6)};
    EXPECT_EQ(expected_product, product);

    std::vector<int> const min = allreduce(send_values, MPI_MIN, mpi);
    std::vector<int> const expected_min = {1, 4};
    EXPECT_EQ(expected_min, min);
}

TEST_F(MPI_BaseLib, AllreduceVectorInplace)
{
    {
        std::vector<int> values = {mpi.rank + 1, mpi.size + mpi.rank + 1};
        allreduceInplace(values, MPI_SUM, mpi);
        std::vector<int> const expected_sum = {(1 + 2 + 3), (4 + 5 + 6)};
        EXPECT_EQ(expected_sum, values);
    }

    {
        std::vector<int> values = {mpi.rank + 1, mpi.size + mpi.rank + 1};
        allreduceInplace(values, MPI_PROD, mpi);
        std::vector<int> const expected_product = {(1 * 2 * 3), (4 * 5 * 6)};
        EXPECT_EQ(expected_product, values);
    }

    {
        std::vector<int> values = {mpi.rank + 1, mpi.size + mpi.rank + 1};
        allreduceInplace(values, MPI_MIN, mpi);
        std::vector<int> const expected_min = {1, 4};
        EXPECT_EQ(expected_min, values);
    }
}

TEST_F(MPI_BaseLib, AnyOf)
{
    EXPECT_FALSE(anyOf(false, mpi));  // false for all ranks.
    EXPECT_TRUE(anyOf(true, mpi));    // true for all ranks.

    // Single rank true and all other false.
    EXPECT_TRUE(anyOf(mpi.rank == 0 ? true : false, mpi));
    // All ranks true and one false.
    EXPECT_TRUE(anyOf(mpi.rank == mpi.size - 1 ? false : true, mpi));
}

TEST_F(MPI_BaseLib, AllOf)
{
    EXPECT_FALSE(allOf(false, mpi));  // false for all ranks.
    EXPECT_TRUE(allOf(true, mpi));    // true for all ranks.

    // Single rank true and all other false.
    EXPECT_FALSE(allOf(mpi.rank == 0 ? true : false, mpi));
    // All ranks true and one false.
    EXPECT_FALSE(allOf(mpi.rank == mpi.size - 1 ? false : true, mpi));
}
#endif

struct ExceptionRuntimeError : std::runtime_error
{
    using std::runtime_error::runtime_error;
};

struct CaptureWarnings
{
    template <typename... Args>
    void operator()(fmt::format_string<Args...> fmt, Args&&... args)
    {
        caught_message = fmt::format(fmt, std::forward<Args>(args)...);
    };

    std::string caught_message;
};

TEST(MPI_BaseLib_Exceptions, NoThrow)
{
    CaptureWarnings warning_callback;
    std::exception_ptr exc;

    EXPECT_NO_THROW(BaseLib::MPI::allRanksThrowOrNone<ExceptionRuntimeError>(
        exc, warning_callback));
    EXPECT_TRUE(warning_callback.caught_message.empty());
}

TEST(MPI_BaseLib_Exceptions, OK)
{
    CaptureWarnings warning_callback;
    std::exception_ptr exc =
        std::make_exception_ptr(ExceptionRuntimeError{"test message"});

    EXPECT_THROW(BaseLib::MPI::allRanksThrowOrNone<ExceptionRuntimeError>(
                     exc, warning_callback),
                 ExceptionRuntimeError);
    EXPECT_TRUE(warning_callback.caught_message.empty());
}

TEST(MPI_BaseLib_Exceptions, OKBase)
{
    CaptureWarnings warning_callback;
    std::exception_ptr exc =
        std::make_exception_ptr(ExceptionRuntimeError{"test message"});

    // Same as above, but with runtime_error as "expected" exception.
    EXPECT_THROW(BaseLib::MPI::allRanksThrowOrNone<std::runtime_error>(
                     exc, warning_callback),
                 ExceptionRuntimeError);
    EXPECT_TRUE(warning_callback.caught_message.empty());
}

TEST(MPI_BaseLib_Exceptions, WarnOtherException)
{
    CaptureWarnings warning_callback;

    std::exception_ptr exc =
        std::make_exception_ptr(ExceptionRuntimeError{"test message"});

    EXPECT_THROW(BaseLib::MPI::allRanksThrowOrNone<std::logic_error>(
                     exc, warning_callback),
                 ExceptionRuntimeError);
    auto const& caught_message = warning_callback.caught_message;
    EXPECT_TRUE(
        caught_message.contains("An exception was thrown on this MPI rank, but "
                                "it's not derived from "))
        << "caught message: " << caught_message;
    EXPECT_TRUE(caught_message.contains(" but rather of type "))
        << "caught message: " << caught_message;
}

TEST(MPI_BaseLib_Exceptions, WarnOtherType)
{
    struct NotAnException
    {
    };

    CaptureWarnings warning_callback;

    std::exception_ptr exc = std::make_exception_ptr(NotAnException{});

    EXPECT_THROW(BaseLib::MPI::allRanksThrowOrNone<std::logic_error>(
                     exc, warning_callback),
                 NotAnException);
    auto const& caught_message = warning_callback.caught_message;
    EXPECT_TRUE(
        caught_message.contains("An exception was thrown on this MPI rank, but "
                                "it's not derived from std::exception."))
        << "caught message: " << caught_message;
}

#ifdef USE_PETSC
TEST(MPI_BaseLib_Exceptions, OnlyOneRankThrows)
{
    struct S : std::runtime_error
    {
        using std::runtime_error::runtime_error;
    };

    Mpi mpi{};
    bool const i_throw = mpi.rank % 3 == 1;

    std::exception_ptr exc =
        i_throw ? std::make_exception_ptr(S{"test message"}) : nullptr;

    bool success = false;
    try
    {
        BaseLib::MPI::allRanksThrowOrNone<std::runtime_error>(exc);
    }
    catch (AnotherMPIRankThrew<std::runtime_error> const& e)
    {
        if (!i_throw)
        {
            success = true;
        }
    }
    catch (std::runtime_error const& e)
    {
        if (i_throw)
        {
            success = true;
        }
    }
    catch (...)
    {
        success = false;
    }

    EXPECT_TRUE(allOf(success));
}
#endif
