/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <limits>
#include <random>
#include <gtest/gtest.h>

#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "ProcessLib/LocalAssemblerInterface.h"
#include "ProcessLib/AnalyticalJacobianAssembler.h"
#include "ProcessLib/CentralDifferencesJacobianAssembler.h"

//! Fills a vector with values whose absolute value is between \c abs_min and
//! \c abs_max.
void fillRandomlyConstrainedAbsoluteValues(std::vector<double>& xs,
                                           double const abs_min,
                                           double const abs_max)
{
    std::random_device rd;
    std::mt19937 random_number_generator(rd());
    double const abs_range = abs_max - abs_min;
    std::uniform_real_distribution<double> rnd(abs_min,
                                               abs_min + 2.0 * abs_range);

    for (auto& x : xs) {
        // v in [ abs_min, abs_min + 2 abs_range ]
        auto v = rnd(random_number_generator);
        if (v > abs_max) {
            // produce negative values
            // (v - abs_range) in [ abs_min, abs_max ]
            v = -(v - abs_range);
        }
        x = v;
    }
}

struct MatDiagX
{
    // M = diag(x1, x2, x3, ...)
    static void getMat(std::vector<double> const& x_data,
                       std::vector<double>& mat_data)
    {
        auto mat =
            MathLib::toZeroedMatrix(mat_data, x_data.size(), x_data.size());
        mat.diagonal().noalias() = MathLib::toVector(x_data, x_data.size());
    }

    // dM/dx * y = diag(y1, y2, y3, ...)
    static void getDMatDxTimesY(std::vector<double> const& x_data,
                                std::vector<double> const& y_data,
                                std::vector<double>& dMdxTy_data)
    {
        auto dMdxTy =
            MathLib::toZeroedMatrix(dMdxTy_data, x_data.size(), x_data.size());
        dMdxTy.diagonal().noalias() = MathLib::toVector(y_data, y_data.size());
    }
};

struct VecX
{
    // v = (x1, x2, x3, ...)
    static void getVec(std::vector<double> const& x_data,
                       std::vector<double>& vec_data)
    {
        vec_data = x_data;
    }

    // dv/dx = Identity
    static void getDVecDx(std::vector<double> const& x_data,
                          std::vector<double>& mat_data)
    {
        auto mat =
            MathLib::toZeroedMatrix(mat_data, x_data.size(), x_data.size());
        mat = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic,
                            Eigen::RowMajor>::Identity(x_data.size(),
                                                       x_data.size());
    }
};

struct MatXY
{
    /*     /  xx  xy  xz ... \
     * M = |  yx  yy  yz ... |
     *     |  zx  zy  zz ... |
     *     \ ... ... ... ... /
     */
    static void getMat(std::vector<double> const& x_data,
                       std::vector<double>& mat_data)
    {
        auto mat =
            MathLib::toZeroedMatrix(mat_data, x_data.size(), x_data.size());
        for (std::size_t r=0; r<x_data.size(); ++r) {
            for (std::size_t c=0; c<x_data.size(); ++c) {
                mat(r, c) = x_data[r] * x_data[c];
            }
        }
    }

    /*             /  xX  xY  xZ ... \
     * dM/dx * y = |  yX  yY  yZ ... | + (x . y) . Identity
     *             |  zX  zY  zZ ... |
     *             \ ... ... ... ... /
     * where x = (x, y, z, ...); y = (X, Y, Z, ...)
     */
    static void getDMatDxTimesY(std::vector<double> const& x_data,
                                std::vector<double> const& y_data,
                                std::vector<double>& dMdxTy_data)
    {
        auto const N = x_data.size();
        auto const x_dot_y =
            MathLib::toVector(x_data, N).dot(MathLib::toVector(y_data, N));
        auto dMdxTy =
            MathLib::toZeroedMatrix(dMdxTy_data, N, N);
        for (std::size_t r=0; r<N; ++r) {
            for (std::size_t c=0; c<N; ++c) {
                dMdxTy(r, c) = x_data[r] * y_data[c];
            }
            dMdxTy(r, r) += x_dot_y;
        }
    }
};

struct VecXRevX
{
    // v = (wz, xy, yx, zw)
    static void getVec(std::vector<double> const& x_data,
                       std::vector<double>& vec_data)
    {
        vec_data = x_data;
        auto const N = x_data.size();
        for (std::size_t i=0; i<N; ++i) {
            vec_data[i] *= x_data[N-i-1];
        }
    }

    /*         / z        w \
     * dv/dx = |    y  x    |
     *         |    y  x    |
     *         \ z        w /
     */
    static void getDVecDx(std::vector<double> const& x_data,
                          std::vector<double>& mat_data)
    {
        auto mat =
            MathLib::toZeroedMatrix(mat_data, x_data.size(), x_data.size());
        auto const N = x_data.size();
        for (std::size_t i=0; i<N; ++i) {
            mat(i, i) += x_data[N-i-1];
            mat(N-i-1, i) += x_data[N-i-1];
        }
    }
};

struct MatDiagXSquared
{
    // M = diag(x^2, y^2, z^2, ...)
    static void getMat(std::vector<double> const& x_data,
                       std::vector<double>& mat_data)
    {
        auto mat =
            MathLib::toZeroedMatrix(mat_data, x_data.size(), x_data.size());

        for (std::size_t i=0; i<x_data.size(); ++i)
            mat(i, i) = x_data[i] * x_data[i];
    }

    // dM/dx * y = diag(2*x*X, 2*y*Y, 2*z*Z, ...)
    static void getDMatDxTimesY(std::vector<double> const& x_data,
                                std::vector<double> const& y_data,
                                std::vector<double>& dMdxTy_data)
    {
        auto dMdxTy =
            MathLib::toZeroedMatrix(dMdxTy_data, x_data.size(), x_data.size());

        for (std::size_t i=0; i<x_data.size(); ++i)
            dMdxTy(i, i) = 2.0 * x_data[i] * y_data[i];
    }
};

struct VecXSquared
{
    // v = (x^2, y^2, z^2, ...)
    static void getVec(std::vector<double> const& x_data,
                       std::vector<double>& vec_data)
    {
        vec_data = x_data;
        for (auto& v : vec_data)
            v *= v;
    }

    // dv/dx = diag(2x, 2y, 2z, ...)
    static void getDVecDx(std::vector<double> const& x_data,
                          std::vector<double>& mat_data)
    {
        auto mat =
            MathLib::toZeroedMatrix(mat_data, x_data.size(), x_data.size());
        for (std::size_t i=0; i< x_data.size(); ++i) {
            mat(i, i) = 2.0 * x_data[i];
        }
    }
};

struct MatXSquaredShifted
{
    /* Tests an asymmetric matrix.
     *
     *     / x^2  y^2  z^2 \
     * M = | z^2  x^2  y^2 |
     *     \ y^2  z^2  x^2 /
     */
    static void getMat(std::vector<double> const& x_data,
                       std::vector<double>& mat_data)
    {
        auto const N = x_data.size();
        auto mat =
            MathLib::toZeroedMatrix(mat_data, N, N);
        for (std::size_t r=0; r<N; ++r) {
            for (std::size_t c=0; c<N; ++c) {
                auto const i = (r + c) % N;
                mat(r, i) = x_data[c] * x_data[c];
            }
        }
    }

    /*             /  2xX  2yY  2zZ \
     * dM/dx * y = |  2xY  2yZ  2zX |
     *             \  2xZ  2yX  2zY /
     * where x = (x, y, z, ...); y = (X, Y, Z, ...)
     */
    static void getDMatDxTimesY(std::vector<double> const& x_data,
                                std::vector<double> const& y_data,
                                std::vector<double>& dMdxTy_data)
    {
        auto const N = x_data.size();
        auto dMdxTy =
            MathLib::toZeroedMatrix(dMdxTy_data, N, N);
        for (std::size_t r=0; r<N; ++r) {
            for (std::size_t c=0; c<N; ++c) {
                auto const i = (r + c) % N;
                dMdxTy(r, c) = 2.0 * x_data[c] * y_data[i];
            }
        }
    }
};

struct MatVecDiagX
{
    using Mat = MatDiagX;
    using Vec = VecX;
    static const double tolerance;
};
const double MatVecDiagX::tolerance = 3e-8;

struct MatVecDiagXSquared
{
    using Mat = MatDiagXSquared;
    using Vec = VecXSquared;
    static const double tolerance;
};
const double MatVecDiagXSquared::tolerance = 2e-7;

struct MatVecXY
{
    using Mat = MatXY;
    using Vec = VecXRevX;
    static const double tolerance;
};
const double MatVecXY::tolerance = 9e-7;

struct MatVecXSquaredShifted
{
    using Mat = MatXSquaredShifted;
    using Vec = VecXRevX;
    static const double tolerance;
};
const double MatVecXSquaredShifted::tolerance = 7e-7;

template <typename MatVec>
class LocalAssemblerM final : public ProcessLib::LocalAssemblerInterface
{
public:
    void assemble(double const /*t*/, std::vector<double> const& local_x,
                  std::vector<double>& local_M_data,
                  std::vector<double>& /*local_K_data*/,
                  std::vector<double>& /*local_b_data*/) override
    {
        MatVec::Mat::getMat(local_x, local_M_data);
    }

    void assembleWithJacobian(double const t,
                              std::vector<double> const& local_x,
                              std::vector<double> const& local_xdot,
                              const double dxdot_dx, const double /*dx_dx*/,
                              std::vector<double>& local_M_data,
                              std::vector<double>& local_K_data,
                              std::vector<double>& local_b_data,
                              std::vector<double>& local_Jac_data) override
    {
        assemble(t, local_x, local_M_data, local_K_data, local_b_data);

        // dM/dx * xdot
        MatVec::Mat::getDMatDxTimesY(local_x, local_xdot, local_Jac_data);

        auto local_Jac =
            MathLib::toMatrix(local_Jac_data, local_x.size(), local_x.size());
        auto local_M =
            MathLib::toMatrix(local_M_data, local_x.size(), local_x.size());
        local_Jac.noalias() += dxdot_dx * local_M;
    }

    static double getTol() { return MatVec::tolerance; }

    static const bool asmM = true;
    static const bool asmK = false;
    static const bool asmb = false;
};

template <typename MatVec>
class LocalAssemblerK final : public ProcessLib::LocalAssemblerInterface
{
public:
    void assemble(double const /*t*/, std::vector<double> const& local_x,
                  std::vector<double>& /*local_M_data*/,
                  std::vector<double>& local_K_data,
                  std::vector<double>& /*local_b_data*/) override
    {
        MatVec::Mat::getMat(local_x, local_K_data);
    }

    void assembleWithJacobian(double const t,
                              std::vector<double> const& local_x,
                              std::vector<double> const& /*local_xdot*/,
                              const double /*dxdot_dx*/, const double dx_dx,
                              std::vector<double>& local_M_data,
                              std::vector<double>& local_K_data,
                              std::vector<double>& local_b_data,
                              std::vector<double>& local_Jac_data) override
    {
        assemble(t, local_x, local_M_data, local_K_data, local_b_data);

        // dK/dx * x
        MatVec::Mat::getDMatDxTimesY(local_x, local_x, local_Jac_data);

        auto local_Jac =
            MathLib::toMatrix(local_Jac_data, local_x.size(), local_x.size());
        auto local_K =
            MathLib::toMatrix(local_K_data, local_x.size(), local_x.size());
        local_Jac.noalias() += dx_dx * local_K;
    }

    static double getTol() { return MatVec::tolerance; }

    static const bool asmM = false;
    static const bool asmK = true;
    static const bool asmb = false;
};

template <typename MatVec>
class LocalAssemblerB final : public ProcessLib::LocalAssemblerInterface
{
public:
    void assemble(double const /*t*/, std::vector<double> const& local_x,
                  std::vector<double>& /*local_M_data*/,
                  std::vector<double>& /*local_K_data*/,
                  std::vector<double>& local_b_data) override
    {
        MatVec::Vec::getVec(local_x, local_b_data);
    }

    void assembleWithJacobian(double const t,
                              std::vector<double> const& local_x,
                              std::vector<double> const& /*local_xdot*/,
                              const double /*dxdot_dx*/, const double /*dx_dx*/,
                              std::vector<double>& local_M_data,
                              std::vector<double>& local_K_data,
                              std::vector<double>& local_b_data,
                              std::vector<double>& local_Jac_data) override
    {
        assemble(t, local_x, local_M_data, local_K_data, local_b_data);

        // db/dx
        MatVec::Vec::getDVecDx(local_x, local_Jac_data);

        auto local_Jac =
            MathLib::toMatrix(local_Jac_data, local_x.size(), local_x.size());
        // J = -db/dx !!
        local_Jac = -local_Jac;
    }

    static double getTol() { return MatVec::tolerance; }

    static const bool asmM = false;
    static const bool asmK = false;
    static const bool asmb = true;
};

template <typename MatVecM, typename MatVecK, typename MatVecB>
class LocalAssemblerMKb final : public ProcessLib::LocalAssemblerInterface
{
public:
    void assemble(double const /*t*/, std::vector<double> const& local_x,
                  std::vector<double>& local_M_data,
                  std::vector<double>& local_K_data,
                  std::vector<double>& local_b_data) override
    {
        MatVecM::Mat::getMat(local_x, local_M_data);
        MatVecK::Mat::getMat(local_x, local_K_data);
        MatVecB::Vec::getVec(local_x, local_b_data);
    }

    void assembleWithJacobian(double const t,
                              std::vector<double> const& local_x,
                              std::vector<double> const& local_xdot,
                              const double dxdot_dx, const double dx_dx,
                              std::vector<double>& local_M_data,
                              std::vector<double>& local_K_data,
                              std::vector<double>& local_b_data,
                              std::vector<double>& local_Jac_data) override
    {
        assemble(t, local_x, local_M_data, local_K_data, local_b_data);

        std::vector<double> local_JacM_data;
        // dM/dx * xdot
        MatVecM::Mat::getDMatDxTimesY(local_x, local_xdot, local_JacM_data);

        std::vector<double> local_JacK_data;
        // dK/dx * x
        MatVecK::Mat::getDMatDxTimesY(local_x, local_x, local_JacK_data);

        std::vector<double> local_JacB_data;
        // db/dx
        MatVecB::Vec::getDVecDx(local_x, local_JacB_data);

        local_Jac_data = local_JacM_data;
        auto const N = local_Jac_data.size();
        for (std::size_t i = 0; i < N; ++i) {
            local_Jac_data[i] += local_JacK_data[i] - local_JacB_data[i];
        }

        auto local_Jac =
            MathLib::toMatrix(local_Jac_data, local_x.size(), local_x.size());
        auto local_M =
            MathLib::toMatrix(local_M_data, local_x.size(), local_x.size());
        auto local_K =
            MathLib::toMatrix(local_K_data, local_x.size(), local_x.size());
        local_Jac.noalias() += dxdot_dx * local_M + dx_dx * local_K;
    }

    static double getTol()
    {
        return std::max(
            {MatVecM::tolerance, MatVecK::tolerance, MatVecB::tolerance});
    }

    static const bool asmM = true;
    static const bool asmK = true;
    static const bool asmb = true;
};

template<class LocAsm>
struct ProcessLibCentralDifferencesJacobianAssembler : public ::testing::Test
{
    static void test()
    {
        // these four local variables will be filled randomly
        std::vector<double> x, xdot;
        double dxdot_dx, dx_dx;

        std::random_device rd;
        std::mt19937 random_number_generator(rd());
        {
            std::uniform_int_distribution<std::size_t> rnd(3, 64);

            auto const size = rnd(random_number_generator);
            x.resize(size);
            xdot.resize(size);

            // all components will be of order of magnitude one
            fillRandomlyConstrainedAbsoluteValues(x, 0.5, 1.5);
            fillRandomlyConstrainedAbsoluteValues(xdot, 0.5, 1.5);
        }
        {
            std::uniform_real_distribution<double> rnd;
            dxdot_dx = rnd(random_number_generator);
            dx_dx = rnd(random_number_generator);
        }

        testInner(x, xdot, dxdot_dx, dx_dx);
    }

private:
    static void testInner(std::vector<double> const& x,
                          std::vector<double> const& xdot,
                          const double dxdot_dx, const double dx_dx)
    {
        ProcessLib::AnalyticalJacobianAssembler jac_asm_ana;
        ProcessLib::CentralDifferencesJacobianAssembler jac_asm_cd({ 1e-8 });
        LocAsm loc_asm;

        double const eps = std::numeric_limits<double>::epsilon();

        std::vector<double> M_data_cd, K_data_cd, b_data_cd, Jac_data_cd,
            M_data_ana, K_data_ana, b_data_ana, Jac_data_ana;
        double const t = 0.0;

        jac_asm_cd.assembleWithJacobian(loc_asm, t, x, xdot, dxdot_dx, dx_dx, M_data_cd,
                                     K_data_cd, b_data_cd, Jac_data_cd);

        jac_asm_ana.assembleWithJacobian(loc_asm, t, x, xdot, dxdot_dx, dx_dx,
                                         M_data_ana, K_data_ana, b_data_ana,
                                         Jac_data_ana);

        if (LocAsm::asmM) {
            ASSERT_EQ(x.size()*x.size(), M_data_cd.size());
            ASSERT_EQ(x.size()*x.size(), M_data_ana.size());
            for (std::size_t i=0; i<x.size()*x.size(); ++i)
                EXPECT_NEAR(M_data_ana[i], M_data_cd[i], eps);
        }

        if (LocAsm::asmK) {
            ASSERT_EQ(x.size()*x.size(), K_data_cd.size());
            ASSERT_EQ(x.size()*x.size(), K_data_ana.size());
            for (std::size_t i=0; i<x.size()*x.size(); ++i)
                EXPECT_NEAR(K_data_ana[i], K_data_cd[i], eps);
        }

        if (LocAsm::asmb) {
            ASSERT_EQ(x.size(), b_data_cd.size());
            ASSERT_EQ(x.size(), b_data_ana.size());
            for (std::size_t i=0; i<x.size(); ++i)
                EXPECT_NEAR(b_data_ana[i], b_data_cd[i], eps);
        }

        ASSERT_EQ(x.size()*x.size(), Jac_data_cd.size());
        ASSERT_EQ(x.size()*x.size(), Jac_data_ana.size());
        for (std::size_t i=0; i<x.size()*x.size(); ++i) {
            // DBUG("%lu, %g, %g", i, Jac_data_ana[i], Jac_data_cd[i]);
            EXPECT_NEAR(Jac_data_ana[i], Jac_data_cd[i], LocAsm::getTol());
        }
    }
};

typedef ::testing::Types<
    // DiagX
    LocalAssemblerM<MatVecDiagX>, LocalAssemblerK<MatVecDiagX>,
    LocalAssemblerB<MatVecDiagX>,
    LocalAssemblerMKb<MatVecDiagX, MatVecDiagX, MatVecDiagX>,
    // DiagXSquared
    LocalAssemblerM<MatVecDiagXSquared>, LocalAssemblerK<MatVecDiagXSquared>,
    LocalAssemblerB<MatVecDiagXSquared>,
    LocalAssemblerMKb<MatVecDiagXSquared, MatVecDiagXSquared,
                      MatVecDiagXSquared>,
    // XY
    LocalAssemblerM<MatVecXY>, LocalAssemblerK<MatVecXY>,
    LocalAssemblerB<MatVecXY>, LocalAssemblerMKb<MatVecXY, MatVecXY, MatVecXY>,
    // XSquaredShifted
    LocalAssemblerM<MatVecXSquaredShifted>,
    LocalAssemblerK<MatVecXSquaredShifted>,
    LocalAssemblerB<MatVecXSquaredShifted>,
    LocalAssemblerMKb<MatVecXSquaredShifted, MatVecXSquaredShifted,
                      MatVecXSquaredShifted>>
    TestCases;

TYPED_TEST_CASE(ProcessLibCentralDifferencesJacobianAssembler, TestCases);

TYPED_TEST(ProcessLibCentralDifferencesJacobianAssembler, Test)
{
    TestFixture::test();
}
