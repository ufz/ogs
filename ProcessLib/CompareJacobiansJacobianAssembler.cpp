/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CompareJacobiansJacobianAssembler.h"

#include <sstream>

#include "CreateJacobianAssembler.h"
#include "MathLib/LinAlg/Eigen/EigenMapTools.h"

namespace
{
//! Dumps a \c double value as a Python script snippet.
void dump_py(std::ostream& fh, std::string const& var, double const val)
{
    fh << var << " = " << val << '\n';
}

//! Dumps a \c std::size_t value as a Python script snippet.
void dump_py(std::ostream& fh, std::string const& var, std::size_t const val)
{
    fh << var << " = " << val << '\n';
}

//! Dumps an arbitrary vector as a Python script snippet.
template <typename Vec>
void dump_py_vec(std::ostream& fh, std::string const& var, Vec const& val)
{
    fh << var << " = np.array([";
    for (decltype(val.size()) i = 0; i < val.size(); ++i)
    {
        if (i != 0)
        {
            if (i % 8 == 0)
            {
                // Print at most eight entries on one line,
                // indent with four spaces.
                fh << ",\n    ";
            }
            else
            {
                fh << ", ";
            }
        }
        fh << val[i];
    }
    fh << "])\n";
}

//! Dumps a \c std::vector<double> as a Python script snippet.
void dump_py(std::ostream& fh, std::string const& var,
             std::vector<double> const& val)
{
    dump_py_vec(fh, var, val);
}

//! Dumps an Eigen vector (array with 1 column) as a Python script snippet.
template <typename Derived>
void dump_py(std::ostream& fh, std::string const& var,
             Eigen::ArrayBase<Derived> const& val,
             std::integral_constant<int, 1> /*unused*/)
{
    dump_py_vec(fh, var, val);
}

//! Dumps an Eigen array as a Python script snippet.
template <typename Derived, int ColsAtCompileTime>
void dump_py(std::ostream& fh, std::string const& var,
             Eigen::ArrayBase<Derived> const& val,
             std::integral_constant<int, ColsAtCompileTime> /*unused*/)
{
    fh << var << " = np.array([\n";
    for (std::ptrdiff_t r = 0; r < val.rows(); ++r)
    {
        if (r != 0)
        {
            fh << ",\n";
        }
        fh << "    [";
        for (std::ptrdiff_t c = 0; c < val.cols(); ++c)
        {
            if (c != 0)
            {
                fh << ", ";
            }
            fh << val(r, c);
        }
        fh << "]";
    }
    fh << "])\n";
}

//! Dumps an Eigen array as a Python script snippet.
template <typename Derived>
void dump_py(std::ostream& fh, std::string const& var,
             Eigen::ArrayBase<Derived> const& val)
{
    dump_py(fh, var, val,
            std::integral_constant<int, Derived::ColsAtCompileTime>{});
}

//! Dumps an Eigen matrix as a Python script snippet.
template <typename Derived>
void dump_py(std::ostream& fh, std::string const& var,
             Eigen::MatrixBase<Derived> const& val)
{
    dump_py(fh, var, val.array());
}

//! Will be printed if some consistency error is detected.
const std::string msg_fatal =
    "The local matrices M or K or the local vectors b assembled with the two "
    "different Jacobian assemblers differ.";

}  // anonymous namespace

namespace ProcessLib
{
void CompareJacobiansJacobianAssembler::assembleWithJacobian(
    LocalAssemblerInterface& local_assembler, double const t, double const dt,
    std::vector<double> const& local_x, std::vector<double> const& local_xdot,
    std::vector<double>& local_M_data, std::vector<double>& local_K_data,
    std::vector<double>& local_b_data, std::vector<double>& local_Jac_data)
{
    ++_counter;

    auto const num_dof = local_x.size();
    auto to_mat = [num_dof](std::vector<double> const& data)
    {
        if (data.empty())
        {
            return Eigen::Map<Eigen::Matrix<
                double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> const>(
                nullptr, 0, 0);
        }

        return MathLib::toMatrix(data, num_dof, num_dof);
    };

    // First assembly -- the one whose results will be added to the global
    // equation system finally.
    _asm1->assembleWithJacobian(local_assembler, t, dt, local_x, local_xdot,
                                local_M_data, local_K_data, local_b_data,
                                local_Jac_data);

    auto const local_M1 = to_mat(local_M_data);
    auto const local_K1 = to_mat(local_K_data);
    auto const local_b1 = MathLib::toVector(local_b_data);

    std::vector<double> local_M_data2;
    std::vector<double> local_K_data2;
    std::vector<double> local_b_data2;
    std::vector<double> local_Jac_data2;

    // Second assembly -- used for checking only.
    _asm2->assembleWithJacobian(local_assembler, t, dt, local_x, local_xdot,
                                local_M_data2, local_K_data2, local_b_data2,
                                local_Jac_data2);

    auto const local_M2 = to_mat(local_M_data2);
    auto const local_K2 = to_mat(local_K_data2);
    auto const local_b2 = MathLib::toVector(local_b_data2);

    auto const local_Jac1 = MathLib::toMatrix(local_Jac_data, num_dof, num_dof);
    auto const local_Jac2 =
        MathLib::toMatrix(local_Jac_data2, num_dof, num_dof);

    auto const abs_diff = (local_Jac2 - local_Jac1).array().eval();
    auto const rel_diff =
        (abs_diff == 0.0)
            .select(abs_diff,
                    2. * abs_diff /
                        (local_Jac1.cwiseAbs() + local_Jac2.cwiseAbs()).array())
            .eval();

    auto const abs_diff_mask =
        (abs_diff.abs() <= _abs_tol)
            .select(decltype(abs_diff)::Zero(abs_diff.rows(), abs_diff.cols()),
                    decltype(abs_diff)::Ones(abs_diff.rows(), abs_diff.cols()))
            .eval();
    auto const rel_diff_mask =
        (rel_diff.abs() <= _rel_tol)
            .select(decltype(rel_diff)::Zero(rel_diff.rows(), rel_diff.cols()),
                    decltype(rel_diff)::Ones(rel_diff.rows(), rel_diff.cols()))
            .eval();

    auto const abs_diff_OK = !abs_diff_mask.any();
    auto const rel_diff_OK = !rel_diff_mask.any();

    std::ostringstream msg_tolerance;
    bool tol_exceeded = true;
    bool fatal_error = false;

    if (abs_diff_OK)
    {
        tol_exceeded = false;
    }
    else
    {
        msg_tolerance << "absolute tolerance of " << _abs_tol << " exceeded";
    }

    if (rel_diff_OK)
    {
        tol_exceeded = false;
    }
    else
    {
        if (!msg_tolerance.str().empty())
        {
            msg_tolerance << " and ";
        }

        msg_tolerance << "relative tolerance of " << _rel_tol << " exceeded";
    }

    // basic consistency check if something went terribly wrong
    auto check_equality = [&fatal_error](auto mat_or_vec1, auto mat_or_vec2)
    {
        if (mat_or_vec1.size() == 0 || mat_or_vec2.size() == 0)
        {
            return;
        }
        if (mat_or_vec1.rows() != mat_or_vec2.rows() ||
            mat_or_vec1.cols() != mat_or_vec2.cols())
        {
            fatal_error = true;
        }
        else if (((mat_or_vec1 - mat_or_vec2).array().cwiseAbs() >
                  std::numeric_limits<double>::epsilon())
                     .any())
        {
            fatal_error = true;
        }
    };

    check_equality(local_M1, local_M2);
    check_equality(local_K1, local_K2);
    check_equality(local_b1, local_b2);

    Eigen::VectorXd res1 = Eigen::VectorXd::Zero(num_dof);
    if (local_M1.size() != 0)
    {
        res1.noalias() += local_M1 * MathLib::toVector(local_xdot);
    }
    if (local_K1.size() != 0)
    {
        res1.noalias() += local_K1 * MathLib::toVector(local_x);
    }
    if (local_b1.size() != 0)
    {
        res1.noalias() -= local_b1;
    }

    Eigen::VectorXd res2 = Eigen::VectorXd::Zero(num_dof);
    if (local_M2.size() != 0)
    {
        res2.noalias() += local_M2 * MathLib::toVector(local_xdot);
    }
    if (local_K2.size() != 0)
    {
        res2.noalias() += local_K2 * MathLib::toVector(local_x);
    }
    if (local_b2.size() != 0)
    {
        res2.noalias() -= local_b2;
    }

    check_equality(res1, res2);

    if (tol_exceeded)
    {
        WARN("Compare Jacobians: {:s}", msg_tolerance.str());
    }

    bool const output = tol_exceeded || fatal_error;

    if (output)
    {
        _log_file << "\n### counter: " << std::to_string(_counter)
                  << " (begin)\n";
    }

    if (fatal_error)
    {
        _log_file << '\n'
                  << "#######################################################\n"
                  << "# FATAL ERROR: " << msg_fatal << '\n'
                  << "#              You cannot expect any meaningful insights "
                     "from the Jacobian data printed below!\n"
                  << "#              The reason for the mentioned differences "
                     "might be\n"
                  << "#              (a) that the assembly routine has side "
                     "effects or\n"
                  << "#              (b) that the assembly routines for M, K "
                     "and b themselves differ.\n"
                  << "#######################################################\n"
                  << '\n';
    }

    if (tol_exceeded)
    {
        _log_file << "# " << msg_tolerance.str() << "\n\n";
    }

    if (output)
    {
        dump_py(_log_file, "counter", _counter);
        dump_py(_log_file, "num_dof", num_dof);
        dump_py(_log_file, "abs_tol", _abs_tol);
        dump_py(_log_file, "rel_tol", _rel_tol);

        _log_file << '\n';

        dump_py(_log_file, "local_x", local_x);
        dump_py(_log_file, "local_x_dot", local_xdot);
        dump_py(_log_file, "dt", dt);

        _log_file << '\n';

        dump_py(_log_file, "Jacobian_1", local_Jac1);
        dump_py(_log_file, "Jacobian_2", local_Jac2);

        _log_file << '\n';

        _log_file << "# Jacobian_2 - Jacobian_1\n";
        dump_py(_log_file, "abs_diff", abs_diff);
        _log_file << "# Componentwise: 2 * abs_diff / (|Jacobian_1| + "
                     "|Jacobian_2|)\n";
        dump_py(_log_file, "rel_diff", rel_diff);

        _log_file << '\n';

        _log_file << "# Masks: 0 ... tolerance met, 1 ... tolerance exceeded\n";
        dump_py(_log_file, "abs_diff_mask", abs_diff_mask);
        dump_py(_log_file, "rel_diff_mask", rel_diff_mask);

        _log_file << '\n';

        dump_py(_log_file, "M_1", local_M1);
        dump_py(_log_file, "M_2", local_M2);
        if (fatal_error && local_M1.size() == local_M2.size())
        {
            dump_py(_log_file, "delta_M", local_M2 - local_M1);
            _log_file << '\n';
        }

        dump_py(_log_file, "K_1", local_K1);
        dump_py(_log_file, "K_2", local_K2);
        if (fatal_error && local_K1.size() == local_K2.size())
        {
            dump_py(_log_file, "delta_K", local_K2 - local_K1);
            _log_file << '\n';
        }

        dump_py(_log_file, "b_1", local_b_data);
        dump_py(_log_file, "b_2", local_b_data2);
        if (fatal_error && local_b1.size() == local_b2.size())
        {
            dump_py(_log_file, "delta_b", local_b2 - local_b1);
            _log_file << '\n';
        }

        dump_py(_log_file, "res_1", res1);
        dump_py(_log_file, "res_2", res2);
        if (fatal_error)
        {
            dump_py(_log_file, "delta_res", res2 - res1);
        }

        _log_file << '\n';

        _log_file << "### counter: " << std::to_string(_counter) << " (end)\n";
    }

    if (fatal_error)
    {
        _log_file << std::flush;
        OGS_FATAL("{:s}", msg_fatal);
    }

    if (tol_exceeded && _fail_on_error)
    {
        _log_file << std::flush;
        OGS_FATAL(
            "OGS failed, because the two Jacobian implementations returned "
            "different results.");
    }
}

std::unique_ptr<CompareJacobiansJacobianAssembler>
createCompareJacobiansJacobianAssembler(BaseLib::ConfigTree const& config)
{
    // TODO doc script corner case: Parameter could occur at different
    // locations.
    //! \ogs_file_param{prj__processes__process__jacobian_assembler__type}
    config.checkConfigParameter("type", "CompareJacobians");

    auto asm1 =
        //! \ogs_file_param{prj__processes__process__jacobian_assembler__CompareJacobians__jacobian_assembler}
        createJacobianAssembler(config.getConfigSubtree("jacobian_assembler"));

    auto asm2 = createJacobianAssembler(
        //! \ogs_file_param{prj__processes__process__jacobian_assembler__CompareJacobians__reference_jacobian_assembler}
        config.getConfigSubtree("reference_jacobian_assembler"));

    //! \ogs_file_param{prj__processes__process__jacobian_assembler__CompareJacobians__abs_tol}
    auto const abs_tol = config.getConfigParameter<double>("abs_tol");
    //! \ogs_file_param{prj__processes__process__jacobian_assembler__CompareJacobians__rel_tol}
    auto const rel_tol = config.getConfigParameter<double>("rel_tol");

    //! \ogs_file_param{prj__processes__process__jacobian_assembler__CompareJacobians__fail_on_error}
    auto const fail_on_error = config.getConfigParameter<bool>("fail_on_error");

    //! \ogs_file_param{prj__processes__process__jacobian_assembler__CompareJacobians__log_file}
    auto const log_file = config.getConfigParameter<std::string>("log_file");

    return std::make_unique<CompareJacobiansJacobianAssembler>(
        std::move(asm1), std::move(asm2), abs_tol, rel_tol, fail_on_error,
        log_file);
}

}  // namespace ProcessLib
