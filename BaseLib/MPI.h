// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <algorithm>
#include <concepts>

#include "Algorithm.h"
#include "DemangleTypeInfo.h"
#include "Error.h"

#ifdef USE_PETSC
#include <mpi.h>
#endif

namespace BaseLib::MPI
{

#ifdef USE_PETSC
extern MPI_Comm OGS_COMM_WORLD;
#endif

struct Setup
{
    Setup(int argc, char* argv[])
    {
#ifdef USE_PETSC
        MPI_Init(&argc, &argv);
#else
        (void)argc;
        (void)argv;
#endif  // USE_PETSC
    }

    ~Setup()
    {
#ifdef USE_PETSC
        MPI_Finalize();
#endif  // USE_PETSC
    }
};

#ifdef USE_PETSC
struct Mpi
{
    Mpi(MPI_Comm const communicator = OGS_COMM_WORLD)
        : communicator(communicator)
    {
        int mpi_init;
        MPI_Initialized(&mpi_init);
        if (mpi_init != 1)
        {
            OGS_FATAL("MPI is not initialized.");
        }
        MPI_Comm_size(communicator, &size);
        MPI_Comm_rank(communicator, &rank);
    }

    MPI_Comm communicator;
    int size;
    int rank;
};

template <typename T>
constexpr MPI_Datatype mpiType()
{
    using U = std::remove_const_t<T>;
    if constexpr (std::is_same_v<U, bool>)
    {
        return MPI_C_BOOL;
    }
    if constexpr (std::is_same_v<U, char>)
    {
        return MPI_CHAR;
    }
    if constexpr (std::is_same_v<U, double>)
    {
        return MPI_DOUBLE;
    }
    if constexpr (std::is_same_v<U, float>)
    {
        return MPI_FLOAT;
    }
    if constexpr (std::is_same_v<U, int>)
    {
        return MPI_INT;
    }
    if constexpr (std::is_same_v<U, std::size_t>)
    {
        return MPI_UNSIGNED_LONG;
    }
    if constexpr (std::is_same_v<U, unsigned int>)
    {
        return MPI_UNSIGNED;
    }
}

template <typename T>
static std::vector<T> allgather(T const& value, Mpi const& mpi)
{
    std::vector<T> result(mpi.size);

    MPI_Allgather(&value, 1, mpiType<T>(), result.data(), 1, mpiType<T>(),
                  mpi.communicator);

    return result;
}

template <typename T>
static std::vector<T> allgather(std::vector<T> const& vector, Mpi const& mpi)
{
    std::size_t const size = vector.size();
    // Flat in memory over all ranks;
    std::vector<T> result(mpi.size * size);

    MPI_Allgather(vector.data(), size, mpiType<T>(), result.data(), size,
                  mpiType<T>(), mpi.communicator);

    return result;
}

template <typename T>
static T allreduce(T const& value, MPI_Op const& mpi_op, Mpi const& mpi)
{
    T result{};

    MPI_Allreduce(&value, &result, 1, mpiType<T>(), mpi_op, mpi.communicator);
    return result;
}

template <typename T>
static std::vector<T> allreduce(std::vector<T> const& vector,
                                MPI_Op const& mpi_op, Mpi const& mpi)
{
    std::size_t const size = vector.size();
    std::vector<T> result(vector.size());

    MPI_Allreduce(vector.data(), result.data(), size, mpiType<T>(), mpi_op,
                  mpi.communicator);
    return result;
}

template <typename T>
static void allreduceInplace(std::vector<T>& vector,
                             MPI_Op const& mpi_op,
                             Mpi const& mpi)
{
    MPI_Allreduce(MPI_IN_PLACE,
                  vector.data(),
                  vector.size(),
                  mpiType<T>(),
                  mpi_op,
                  mpi.communicator);
}

/// Allgather for variable data. Returns offsets in the receive buffer.
/// The receive buffer is resized to accommodate the gathered data.
template <typename T>
static std::vector<int> allgatherv(
    std::span<T> const send_buffer,
    std::vector<std::remove_const_t<T>>& receive_buffer,
    Mpi const& mpi)
{
    // Determine the number of elements to send
    int const size = static_cast<int>(send_buffer.size());

    // Gather sizes from all ranks
    std::vector<int> const sizes = allgather(size, mpi);

    // Compute offsets based on counts
    std::vector<int> const offsets = BaseLib::sizesToOffsets(sizes);

    // Resize receive buffer to hold all gathered data
    receive_buffer.resize(offsets.back());

    MPI_Allgatherv(send_buffer.data(), size, mpiType<T>(),
                   receive_buffer.data(), sizes.data(), offsets.data(),
                   mpiType<T>(), mpi.communicator);

    return offsets;
}
#endif

/// The reduction is implemented transparently for with and without MPI. In the
/// latter case the input value is returned.
static inline bool anyOf(bool const val
#ifdef USE_PETSC
                         ,
                         Mpi const& mpi = Mpi{OGS_COMM_WORLD}
#endif
)
{
#ifdef USE_PETSC
    return allreduce(val, MPI_LOR, mpi);
#else
    return val;
#endif
}

/// The reduction is implemented transparently for with and without MPI. In the
/// latter case the input value is returned.
static inline bool allOf(bool const val
#ifdef USE_PETSC
                         ,
                         Mpi const& mpi = Mpi{OGS_COMM_WORLD}
#endif
)
{
    return !anyOf(!val
#ifdef USE_PETSC
                  ,
                  mpi
#endif
    );
}

/// Class indicating that another MPI rank threw an exception.
///
/// This class is derived from BaseException. Hence, it can be caught with \code
/// catch(BaseException const&) \endcode. That makes it possible to handle both
/// \code BaseException \endcode and \code AnotherMPIRankThrew<BaseException>
/// \endcode in the same way.
template <typename BaseException>
    requires std::derived_from<BaseException, std::exception> &&
             (  // The used ctor excludes std::exception itself
                 !std::same_as<BaseException, std::exception>)
class AnotherMPIRankThrew : public BaseException
{
public:
    using BaseException::BaseException;

    AnotherMPIRankThrew()
        : BaseException{"Another MPI rank threw an exception."}
    {
    }
};

/// This function throws if and only if the passed exception is non-null on any
/// MPI rank.
///
/// The thrown exception is:
/// - the passed one on any MPI rank where the passed exception is non-null
/// - AnotherMPIRankThrew<Exception> on all other MPI ranks
///
/// This function also checks if the thrown exception is derived from \c
/// Exception.
template <typename Exception>
void allRanksThrowOrNone([[maybe_unused]] std::exception_ptr const& exception,
                         auto&& warning_callback)
{
    // std::exception would lead to duplicate catch clauses below and would not
    // work together with the current implementation of AnotherMPIRankThrew.
    static_assert(!std::is_same_v<Exception, std::exception>);

    bool const exception_was_thrown = anyOf(exception != nullptr);

    [[unlikely]] if (exception_was_thrown)
    {
        if (exception)
        {
            try
            {
                std::rethrow_exception(exception);
            }
            catch (Exception const&)
            {
                // OK. Argument exception is derived from class Exception.
                throw;
            }
            catch (std::exception const& e)
            {
                warning_callback(
                    "An exception was thrown on this MPI rank, but it's not "
                    "derived from {}, but rather of type {}",
                    BaseLib::typeToString<Exception>(),
                    BaseLib::demangle(
                        typeid(e).name() /* demangle the runtime type of e */));
                throw;
            }
            catch (...)
            {
                warning_callback(
                    "An exception was thrown on this MPI rank, but it's not "
                    "derived from std::exception.");
                throw;
            }
        }

        throw AnotherMPIRankThrew<Exception>{};
    }
}

/// Same as above but using the standard WARN() function for warnings.
template <typename Exception>
void allRanksThrowOrNone([[maybe_unused]] std::exception_ptr const& exception)
{
    auto warning_callback =
        []<typename... Args>(fmt::format_string<Args...> fmt, Args&&... args)
    { WARN(fmt, std::forward<Args>(args)...); };

    allRanksThrowOrNone<Exception>(exception, warning_callback);
}

}  // namespace BaseLib::MPI
