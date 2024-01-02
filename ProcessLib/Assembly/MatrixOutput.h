/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <fstream>
#include <functional>
#include <mutex>
#include <unordered_set>

#include "MathLib/LinAlg/GlobalMatrixVectorTypes.h"

namespace ProcessLib::Assembly
{
//! Writes global matrices to disk for debugging purposes.
struct GlobalMatrixOutput
{
    GlobalMatrixOutput();

    void operator()(double const t, int const process_id, GlobalMatrix const& M,
                    GlobalMatrix const& K, GlobalVector const& b,
                    GlobalMatrix const* const Jac = nullptr);

private:
    std::string filenamePrefix_;

    //! Used to distinguish output files of global matrices at the same time,
    //! but, e.g., at different non-linear iterations etc.
    std::size_t counter_ = 0;

    bool do_output_ = false;
};

//! Writes local matrices to disk for debugging purposes.
struct LocalMatrixOutput
{
    LocalMatrixOutput();

    void operator()(double const t, int const process_id,
                    std::size_t const element_id,
                    std::vector<double> const& local_M_data,
                    std::vector<double> const& local_K_data,
                    std::vector<double> const& local_b_data,
                    std::vector<double> const* const local_Jac_data = nullptr);

    ~LocalMatrixOutput();

private:
    bool isOutputRequested(std::size_t const element_id) const;

    std::mutex mutex_;
    std::ofstream outputFile_;
    std::function<bool(std::size_t)> output_element_predicate_;
};

namespace detail
{
// helper function, exposed for unit testing, only
std::function<bool(std::size_t)> createLocalMatrixOutputElementPredicate(
    std::string const& element_ids_str);

// helper function, exposed for unit testing, only
std::unordered_set<std::size_t> parseSetOfSizeT(std::string const& str,
                                                std::string const& warn_msg);
}  // namespace detail
}  // namespace ProcessLib::Assembly
