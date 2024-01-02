/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "MatrixOutput.h"

#include <spdlog/fmt/bundled/ostream.h>

#include <optional>
#include <unordered_set>

#include "BaseLib/Error.h"
#include "BaseLib/Logging.h"
#include "MathLib/FormattingUtils.h"
#include "MathLib/LinAlg/Eigen/EigenMapTools.h"

namespace
{
std::string getSeparatorAfterFilenamePrefix(std::string const& filenamePrefix)
{
    return filenamePrefix.empty() || filenamePrefix.ends_with('/') ||
                   filenamePrefix.ends_with('\\')
               ? ""
               : "_";
}

#ifndef USE_PETSC
static void outputGlobalMatrix(GlobalMatrix const& mat, std::ostream& os)
{
    os << std::setprecision(16) << "(" << mat.getNumberOfRows() << " x "
       << mat.getNumberOfColumns() << ")\n";
    mat.write(os);
}

static void outputGlobalVector(GlobalVector const& vec, std::ostream& os)
{
    os << std::setprecision(16) << "(" << vec.size() << ")\n";
    os << vec.getRawVector() << '\n';
}

std::ofstream openGlobalMatrixOutputFile(std::string const& filenamePrefix,
                                         std::size_t const counter,
                                         double const t, int const process_id,
                                         std::string const& which_matrix,
                                         std::string const& extension)
{
    auto const filename = fmt::format(
        "{}{}ogs_global_matrix_cnt_{:03}_t_{:g}_pcs_{}_{}.{}", filenamePrefix,
        getSeparatorAfterFilenamePrefix(filenamePrefix), counter, t, process_id,
        which_matrix, extension);

    std::ofstream fh{filename};

    if (!fh)
    {
        OGS_FATAL("Could not open file `{}' for global matrix debug output",
                  filename);
    }

    return fh;
}
#endif

static std::optional<std::string> getEnvironmentVariable(
    std::string const& env_var)
{
    char const* const prefix = std::getenv(env_var.c_str());
    return prefix ? std::make_optional(prefix) : std::nullopt;
}

std::string localMatrixOutputFilename(std::string const& filenamePrefix)
{
    return fmt::format("{}{}ogs_local_matrix.log", filenamePrefix,
                       getSeparatorAfterFilenamePrefix(filenamePrefix));
}

Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic,
                               Eigen::RowMajor>>
toSquareMatrixRowMajor(std::vector<double> entries)
{
    auto const num_r_c =
        static_cast<Eigen::Index>(std::round(std::sqrt(entries.size())));

    return {entries.data(), num_r_c, num_r_c};
}
}  // namespace

namespace ProcessLib::Assembly
{
namespace detail
{
std::unordered_set<std::size_t> parseSetOfSizeT(std::string const& str,
                                                std::string const& warn_msg)
{
    std::istringstream sstr{str};
    std::unordered_set<std::size_t> result;
    std::ptrdiff_t value;

    while (sstr >> value)
    {
        [[likely]] if (value >= 0)
        {
            result.insert(value);
        }
        else
        {
            if (!warn_msg.empty())
            {
                WARN("{}", warn_msg);
            }

            return result;
        }
    }

    if (!sstr.eof() && !warn_msg.empty())  // The stream is not read until
                                           // the end, must be an error.
    {
        WARN("{}", warn_msg);
    }

    return result;
}

std::function<bool(std::size_t)> createLocalMatrixOutputElementPredicate(
    std::string const& element_ids_str)
{
    if (element_ids_str.empty())
    {
        return {};
    }

    if (element_ids_str == "*")
    {
        return [](std::size_t) { return true; };
    }

    auto element_ids = parseSetOfSizeT(
        element_ids_str,
        "Error parsing list of element ids for local matrix debug "
        "output. We'll try to proceed anyway, as best as we can.");

    if (element_ids.empty())
    {
        return {};
    }

    return [element_ids = std::move(element_ids)](std::size_t element_id)
    { return element_ids.contains(element_id); };
}
}  // namespace detail

GlobalMatrixOutput::GlobalMatrixOutput()
{
    auto opt_prefix = getEnvironmentVariable("OGS_GLOBAL_MAT_OUT_PREFIX");

    if (opt_prefix.has_value())
    {
#ifndef USE_PETSC
        do_output_ = true;
        filenamePrefix_ = std::move(*opt_prefix);
#else
        // TODO implement. PETScMatrix's viewer() method could be used for that.
        WARN(
            "You requested global matrix output (OGS_GLOBAL_MAT_OUT_PREFIX is "
            "set), which is not yet implemented for PETSc matrices.");
#endif
    }
}

void GlobalMatrixOutput::operator()(double const t, int const process_id,
                                    GlobalMatrix const& M,
                                    GlobalMatrix const& K,
                                    GlobalVector const& b,
                                    GlobalMatrix const* const Jac)
{
    if (!do_output_)
    {
        return;
    }

#ifndef USE_PETSC
    ++counter_;

    {
        auto fh = openGlobalMatrixOutputFile(filenamePrefix_, counter_, t,
                                             process_id, "M", "mat");

        fh << "M ";
        outputGlobalMatrix(M, fh);
    }

    {
        auto fh = openGlobalMatrixOutputFile(filenamePrefix_, counter_, t,
                                             process_id, "K", "mat");

        fh << "K ";
        outputGlobalMatrix(K, fh);
    }

    {
        auto fh = openGlobalMatrixOutputFile(filenamePrefix_, counter_, t,
                                             process_id, "b", "vec");

        fh << "b ";
        outputGlobalVector(b, fh);
    }

    if (Jac)
    {
        auto fh = openGlobalMatrixOutputFile(filenamePrefix_, counter_, t,
                                             process_id, "Jac", "mat");

        fh << "Jac ";
        outputGlobalMatrix(*Jac, fh);
    }
#else
    // do nothing, warning message already printed in the constructor
    (void)t;
    (void)process_id;
    (void)M;
    (void)K;
    (void)b;
    (void)Jac;
#endif
}

LocalMatrixOutput::LocalMatrixOutput()
{
    auto const opt_prefix = getEnvironmentVariable("OGS_LOCAL_MAT_OUT_PREFIX");
    auto const opt_elements =
        getEnvironmentVariable("OGS_LOCAL_MAT_OUT_ELEMENTS");

    if (!opt_prefix)
    {
        if (opt_elements)
        {
            WARN(
                "Environment variable OGS_LOCAL_MAT_OUT_ELEMENTS is set, but "
                "OGS_LOCAL_MAT_OUT_PREFIX is not. Local matrix debug output "
                "will be disabled.");
        }

        return;
    }
    if (!opt_elements)
    {
        WARN(
            "Environment variable OGS_LOCAL_MAT_OUT_PREFIX is set, but "
            "OGS_LOCAL_MAT_OUT_ELEMENTS is not. Local matrix debug output "
            "will be disabled.");
        return;
    }

    output_element_predicate_ =
        detail::createLocalMatrixOutputElementPredicate(*opt_elements);

    if (!output_element_predicate_)
    {
        WARN(
            "Environment variable OGS_LOCAL_MAT_OUT_ELEMENTS not set "
            "properly. Local matrix debug output will be disabled.");
        return;
    }

    auto const outputFilename = localMatrixOutputFilename(*opt_prefix);
    outputFile_.open(outputFilename);

    if (!outputFile_)
    {
        OGS_FATAL(
            "File `{}' for local matrix debug output could not be "
            "opened.",
            outputFilename);
    }

    DBUG("Successfully opened local matrix debug output file {}.",
         outputFilename);
}

void LocalMatrixOutput::operator()(
    double const t, int const process_id, std::size_t const element_id,
    std::vector<double> const& local_M_data,
    std::vector<double> const& local_K_data,
    std::vector<double> const& local_b_data,
    std::vector<double> const* const local_Jac_data)
{
    [[likely]] if (!isOutputRequested(element_id))
    {
        return;
    }

    std::lock_guard lock_guard{mutex_};

    auto& fh = outputFile_;

    DBUG("Writing to local matrix debug output file...");

    fmt::print(fh, "## t = {:.15g}, process id = {}, element id = {}\n\n", t,
               process_id, element_id);

    if (!local_M_data.empty())
    {
        DBUG("... M");
        fmt::print(fh, "# M\n{}\n\n", toSquareMatrixRowMajor(local_M_data));
    }

    if (!local_K_data.empty())
    {
        DBUG("... K");
        fmt::print(fh, "# K\n{}\n\n", toSquareMatrixRowMajor(local_K_data));
    }

    if (!local_b_data.empty())
    {
        DBUG("... b");
        fmt::print(fh, "# b\n{}\n\n", MathLib::toVector(local_b_data));
    }

    if (local_Jac_data && !local_Jac_data->empty())
    {
        DBUG("... Jac");
        fmt::print(fh, "# Jac\n{}\n\n\n",
                   toSquareMatrixRowMajor(*local_Jac_data));
    }
}

bool LocalMatrixOutput::isOutputRequested(std::size_t const element_id) const
{
    return output_element_predicate_ && output_element_predicate_(element_id);
}

LocalMatrixOutput::~LocalMatrixOutput() = default;
}  // namespace ProcessLib::Assembly
