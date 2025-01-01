/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <exception>
#include <mutex>

/// Initial code from
/// https://stackoverflow.com/questions/11828539/elegant-exception-handling-in-openmp
///
/// The implementation does not guarantee, that while the exception is being
/// checked for nullptr in the bool-conversion operator, the capture is not
/// changing the pointer. This might lead to an exception being lost if there
/// are two exceptions thrown simultaneously.
class ThreadException
{
public:
    void capture()
    {
        std::unique_lock<std::mutex> guard{lock_};
        exception_ = std::current_exception();
    }

    void rethrow()
    {
        if (exception_)
        {
            std::rethrow_exception(exception_);
        }
    }

    explicit operator bool() const noexcept { return exception_ != nullptr; }

private:
    std::exception_ptr exception_ = nullptr;
    std::mutex lock_;
};
