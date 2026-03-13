// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "CompareJacobiansJacobianAssembler.h"

#include <fstream>
#include <limits>
#include <sstream>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "BaseLib/ConfigTree.h"
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

//! Absolute and (symmetric) relative difference as Eigen::Array
auto absDiffRelDiff(auto const& mat_or_vec1, auto const& mat_or_vec2)
{
    auto abs_diff = (mat_or_vec2 - mat_or_vec1).array().eval();
    auto const rel_diff =
        (abs_diff == 0.0)
            .select(
                abs_diff,
                2. * abs_diff /
                    (mat_or_vec1.cwiseAbs() + mat_or_vec2.cwiseAbs()).array())
            .eval();
    return std::pair{std::move(abs_diff), std::move(rel_diff)};
}

template <typename MatOrVec>
auto oneIfAboveThresholdElseZero(MatOrVec&& mat_or_vec, double const threshold)
{
    return (mat_or_vec <= threshold)
        .select(MatOrVec::Zero(mat_or_vec.rows(), mat_or_vec.cols()),
                MatOrVec::Ones(mat_or_vec.rows(), mat_or_vec.cols()))
        .eval();
}

// basic consistency check if something went terribly wrong
bool isSimilar(auto const& mat_or_vec1, auto const& mat_or_vec2,
               double const abs_tol, double const rel_tol)
{
    if (mat_or_vec1.size() == 0 || mat_or_vec2.size() == 0)
    {
        // either size is 0, ignore
        return true;
    }
    if (mat_or_vec1.rows() != mat_or_vec2.rows() ||
        mat_or_vec1.cols() != mat_or_vec2.cols())
    {
        return false;
    }

    auto const [abs_diff, rel_diff] = absDiffRelDiff(mat_or_vec1, mat_or_vec2);
    auto const abs_tol_exceeded = abs_diff > abs_tol;
    auto const rel_tol_exceeded = rel_diff > rel_tol;

    // similar if for no entry abs and rel tols are exceeded at the same time
    return !(abs_tol_exceeded && rel_tol_exceeded).any();
}

}  // anonymous namespace

namespace ProcessLib
{
namespace detail
{
//! Data shared among all copied CompareJacobiansJacobianAssembler instances
struct CompareJacobiansJacobianAssemblerImpl
{
    friend class ProcessLib::CompareJacobiansJacobianAssembler;

    CompareJacobiansJacobianAssemblerImpl(
        std::unique_ptr<AbstractJacobianAssembler>&& asm1,
        std::unique_ptr<AbstractJacobianAssembler>&& asm2,
        double abs_tol_Jac,
        double rel_tol_Jac,
        double abs_tol_res,
        double rel_tol_res,
        bool fail_on_error,
        std::string const& log_file_path)
        : _asm1{std::move(asm1)},
          _asm2{std::move(asm2)},
          _abs_tol_Jac{abs_tol_Jac},
          _rel_tol_Jac{rel_tol_Jac},
          _abs_tol_res{abs_tol_res},
          _rel_tol_res{rel_tol_res},
          _fail_on_error{fail_on_error},
          _log_file{log_file_path}
    {
        _log_file.precision(std::numeric_limits<double>::max_digits10);
        _log_file << "#!/usr/bin/env python\n"
                     "import numpy as np\n"
                     "from numpy import nan\n"
                  << std::endl;
    }

    void assembleWithJacobian(LocalAssemblerInterface& local_assembler,
                              double const t, double const dt,
                              std::vector<double> const& local_x,
                              std::vector<double> const& local_x_prev,
                              std::vector<double>& local_b_data,
                              std::vector<double>& local_Jac_data);

private:
    std::unique_ptr<AbstractJacobianAssembler> _asm1;
    std::unique_ptr<AbstractJacobianAssembler> _asm2;

    double const _abs_tol_Jac;
    double const _rel_tol_Jac;
    double const _abs_tol_res;
    double const _rel_tol_res;

    //! Whether to abort if the tolerances are exceeded.
    bool const _fail_on_error;

    //! Path where a Python script will be placed, which contains
    //! information about exceeded tolerances and assembled local
    //! matrices.
    std::ofstream _log_file;

    //! Counter used for identifying blocks in the \c _log_file. It is
    //! incremented upon each call of the assembly routine, i.e., for
    //! each element in each iteration etc.
    std::size_t _counter = 0;
};

void CompareJacobiansJacobianAssemblerImpl::assembleWithJacobian(
    LocalAssemblerInterface& local_assembler, double const t, double const dt,
    std::vector<double> const& local_x, std::vector<double> const& local_x_prev,
    std::vector<double>& local_b_data, std::vector<double>& local_Jac_data)
{
    ++_counter;

    auto const num_dof = local_x.size();

    // First assembly -- the one whose results will be added to the global
    // equation system finally.
    _asm1->assembleWithJacobian(local_assembler, t, dt, local_x, local_x_prev,
                                local_b_data, local_Jac_data);

    auto const local_b1 = MathLib::toVector(local_b_data);

    std::vector<double> local_b_data2;
    std::vector<double> local_Jac_data2;

    // Second assembly -- used for checking only.
    _asm2->assembleWithJacobian(local_assembler, t, dt, local_x, local_x_prev,
                                local_b_data2, local_Jac_data2);

    auto const local_b2 = MathLib::toVector(local_b_data2);

    auto const local_Jac1 = MathLib::toMatrix(local_Jac_data, num_dof, num_dof);
    auto const local_Jac2 =
        MathLib::toMatrix(local_Jac_data2, num_dof, num_dof);

    auto const [abs_diff, rel_diff] = absDiffRelDiff(local_Jac1, local_Jac2);

    auto const abs_diff_mask =
        oneIfAboveThresholdElseZero(abs_diff.abs(), _abs_tol_Jac);
    auto const rel_diff_mask =
        oneIfAboveThresholdElseZero(rel_diff.abs(), _rel_tol_Jac);

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
        msg_tolerance << "absolute tolerance of " << _abs_tol_Jac
                      << " exceeded";
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

        msg_tolerance << "relative tolerance of " << _rel_tol_Jac
                      << " exceeded";
    }

    fatal_error |= !isSimilar(local_b1, local_b2, _abs_tol_res, _rel_tol_res);

    Eigen::VectorXd res1 = Eigen::VectorXd::Zero(num_dof);
    auto const x = MathLib::toVector(local_x);
    auto const x_dot = ((x - MathLib::toVector(local_x_prev)) / dt).eval();
    if (local_b1.size() != 0)
    {
        res1.noalias() -= local_b1;
    }

    Eigen::VectorXd res2 = Eigen::VectorXd::Zero(num_dof);
    if (local_b2.size() != 0)
    {
        res2.noalias() -= local_b2;
    }

    fatal_error |= !isSimilar(res1, res2, _abs_tol_res, _rel_tol_res);

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
                  << "#              (b) that the assembly routines for b "
                     "themselves differ.\n"
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
        dump_py(_log_file, "abs_tol", _abs_tol_Jac);
        dump_py(_log_file, "rel_tol", _rel_tol_Jac);

        _log_file << '\n';

        dump_py(_log_file, "local_x", local_x);
        dump_py(_log_file, "local_x_prev", local_x_prev);
        dump_py(_log_file, "dt", dt);

        _log_file << '\n';

        dump_py(_log_file, "Jacobian_1", local_Jac1);
        dump_py(_log_file, "Jacobian_2", local_Jac2);

        _log_file << '\n';

        _log_file << "# Jacobian_2 - Jacobian_1\n";
        dump_py(_log_file, "abs_diff", abs_diff);
        _log_file << "# max(|abs_diff|) = " << abs_diff.abs().maxCoeff()
                  << '\n';
        _log_file << "# Componentwise: 2 * abs_diff / (|Jacobian_1| + "
                     "|Jacobian_2|)\n";
        dump_py(_log_file, "rel_diff", rel_diff);
        _log_file << "# max(|rel_diff|) = " << rel_diff.abs().maxCoeff()
                  << '\n';

        _log_file << '\n';

        _log_file << "# Masks: 0 ... tolerance met, 1 ... tolerance exceeded\n";
        dump_py(_log_file, "abs_diff_mask", abs_diff_mask);
        dump_py(_log_file, "rel_diff_mask", rel_diff_mask);

        _log_file << '\n';

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
}  // namespace detail

CompareJacobiansJacobianAssembler::CompareJacobiansJacobianAssembler(
    std::unique_ptr<AbstractJacobianAssembler>&& asm1,
    std::unique_ptr<AbstractJacobianAssembler>&& asm2, double abs_tol_Jac,
    double rel_tol_Jac, double abs_tol_res, double rel_tol_res,
    bool fail_on_error, std::string const& log_file_path)
    : impl_{std::make_shared<detail::CompareJacobiansJacobianAssemblerImpl>(
          std::move(asm1), std::move(asm2), abs_tol_Jac, rel_tol_Jac,
          abs_tol_res, rel_tol_res, fail_on_error, log_file_path)}
{
}

CompareJacobiansJacobianAssembler::CompareJacobiansJacobianAssembler(
    std::shared_ptr<detail::CompareJacobiansJacobianAssemblerImpl> impl,
    CompareJacobiansJacobianAssembler::Key)
    : impl_{std::move(impl)}
{
}

void CompareJacobiansJacobianAssembler::assembleWithJacobian(
    LocalAssemblerInterface& local_assembler, double const t, double const dt,
    std::vector<double> const& local_x, std::vector<double> const& local_x_prev,
    std::vector<double>& local_b_data, std::vector<double>& local_Jac_data)
{
    impl_->assembleWithJacobian(local_assembler, t, dt, local_x, local_x_prev,
                                local_b_data, local_Jac_data);
}

std::unique_ptr<AbstractJacobianAssembler>
CompareJacobiansJacobianAssembler::copy() const
{
#ifdef _OPENMP
    if (omp_get_thread_num() != 0)
    {
        OGS_FATAL(
            "CompareJacobiansJacobianAssembler cannot be used concurrently. "
            "Please restrict yourself to one assembly thread "
            "(OGS_ASM_THREADS=1).");
    }
#endif

    return std::make_unique<CompareJacobiansJacobianAssembler>(impl_, Key{});
}

void CompareJacobiansJacobianAssembler::checkPerturbationSize(
    int const max_non_deformation_dofs_per_node) const
{
    impl_->_asm1->checkPerturbationSize(max_non_deformation_dofs_per_node);
    impl_->_asm2->checkPerturbationSize(max_non_deformation_dofs_per_node);
}

void CompareJacobiansJacobianAssembler::setNonDeformationComponentIDs(
    std::vector<int> const& non_deformation_component_ids)
{
    impl_->_asm1->setNonDeformationComponentIDs(non_deformation_component_ids);
    impl_->_asm2->setNonDeformationComponentIDs(non_deformation_component_ids);
}

void CompareJacobiansJacobianAssembler::
    setNonDeformationComponentIDsNoSizeCheck(
        std::vector<int> const& non_deformation_component_ids)
{
    impl_->_asm1->setNonDeformationComponentIDsNoSizeCheck(
        non_deformation_component_ids);
    impl_->_asm2->setNonDeformationComponentIDsNoSizeCheck(
        non_deformation_component_ids);
}

bool CompareJacobiansJacobianAssembler::needsPicardAssembly() const
{
    return impl_->_asm1->needsPicardAssembly() ||
           impl_->_asm2->needsPicardAssembly();
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

    //! \ogs_file_param{prj__processes__process__jacobian_assembler__CompareJacobians__abs_tol_res}
    auto const abs_tol_res = config.getConfigParameter<double>("abs_tol_res");
    //! \ogs_file_param{prj__processes__process__jacobian_assembler__CompareJacobians__rel_tol_res}
    auto const rel_tol_res = config.getConfigParameter<double>("rel_tol_res");

    //! \ogs_file_param{prj__processes__process__jacobian_assembler__CompareJacobians__fail_on_error}
    auto const fail_on_error = config.getConfigParameter<bool>("fail_on_error");

    //! \ogs_file_param{prj__processes__process__jacobian_assembler__CompareJacobians__log_file}
    auto const log_file = config.getConfigParameter<std::string>("log_file");

    return std::make_unique<CompareJacobiansJacobianAssembler>(
        std::move(asm1), std::move(asm2), abs_tol, rel_tol, abs_tol_res,
        rel_tol_res, fail_on_error, log_file);
}
}  // namespace ProcessLib
