/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * The code of this file is used to decouple the evaluation of matrix elements from the rest of OGS6,
 * not all of OGS6 has to be recompiled every time a small change is done.
 */

#ifndef PROCESS_LIB_TES_FEM_NOTPL_H_
#define PROCESS_LIB_TES_FEM_NOTPL_H_

#include <memory>
#include <Eigen/Eigen>

#include "TESProcess-notpl.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"

#include "VariableTransformation.h"

namespace ProcessLib
{

namespace TES
{

enum class SecondaryVariables {
    SOLID_DENSITY, REACTION_RATE,
    VELOCITY_X, VELOCITY_Y, VELOCITY_Z,
    VAPOUR_PARTIAL_PRESSURE,
    RELATIVE_HUMIDITY,
    LOADING,
    EQUILIBRIUM_LOADING,
    REACTION_DAMPING_FACTOR
};


template<typename ShpPol, unsigned NIntPts, unsigned NodalDOF, unsigned Dim>
struct DataTraitsFixed
{
    using ShapeMatrices = typename ShpPol::ShapeMatrices;

    template<int N, int M>
    using Mat = typename ShpPol::template MatrixType<N, M>;
    template<int N>
    using Vec = typename ShpPol::template VectorType<N>;

    using MatrixDimDim = Mat<Dim, Dim>;

    using LocalMatrix = Mat<NIntPts*NodalDOF, NIntPts*NodalDOF>;
    using LocalVector = Vec<NIntPts*NodalDOF>;

    using Vector1Comp = typename ShapeMatrices::ShapeType;

    using LaplaceMatrix = Mat<Dim*NodalDOF, Dim*NodalDOF>;


    // block dim x dim, fixed-size const matrix
    // TODO: swap variable names width <--> height
    template<typename Mat>
    static
    typename std::enable_if<NodalDOF != 0,
        decltype(std::declval<const Mat>().template block<Dim, Dim>(0u, 0u))
    >::type
    blockDimDim(Mat const& mat,
                unsigned top, unsigned left, unsigned width, unsigned height)
    {
        assert(width==Dim && height==Dim);
        (void) width; (void) height;
        return mat.template block<Dim, Dim>(top, left);
    }
    // block dim x dim, dynamic-size const matrix
    template<typename Mat>
    static
    typename std::enable_if<NodalDOF == 0,
        decltype(std::declval<const Mat>().block(0u, 0u, 0u, 0u))
    >::type
    blockDimDim(Mat const& mat,
                unsigned top, unsigned left, unsigned width, unsigned height)
    {
        assert(width == height);
        return mat.block(top, left, width, height);
    }
    // block dim x dim, fixed-size non-const matrix
    template<typename Mat>
    static
    typename std::enable_if<NodalDOF != 0,
        decltype(std::declval<Mat>().template block<Dim, Dim>(0u, 0u))
    >::type
    blockDimDim(Mat& mat,
                unsigned top, unsigned left, unsigned width, unsigned height)
    {
        assert(width==Dim && height==Dim);
        (void) width; (void) height;
        return mat.template block<Dim, Dim>(top, left);
    }
    // block dim x dim, dynamic-size non-const matrix
    template<typename Mat>
    static
    typename std::enable_if<NodalDOF == 0,
        decltype(std::declval<Mat>().block(0u, 0u, 0u, 0u))
    >::type
    blockDimDim(Mat& mat,
                unsigned top, unsigned left, unsigned width, unsigned height)
    {
        assert(width == height);
        return mat.block(top, left, width, height);
    }

    // block gauss x gauss, fixed-size const matrix
    template<typename Mat>
    static
    typename std::enable_if<NodalDOF != 0,
        decltype(std::declval<const Mat>().template block<NIntPts, NIntPts>(0u, 0u))
    >::type
    blockShpShp(Mat const& mat,
                unsigned top, unsigned left, unsigned width, unsigned height)
    {
        assert(width==NIntPts && height==NIntPts);
        (void) width; (void) height;
        return mat.template block<NIntPts, NIntPts>(top, left);
    }
    // block gauss x gauss, dynamic-size const matrix
    template<typename Mat>
    static
    typename std::enable_if<NodalDOF == 0,
        decltype(std::declval<const Mat>().block(0u, 0u, 0u, 0u))
    >::type
    blockShpShp(Mat const& mat,
                unsigned top, unsigned left, unsigned width, unsigned height)
    {
        assert(width == height);
        return mat.block(top, left, width, height);
    }
    // block gauss x gauss, fixed-size non-const matrix
    template<typename Mat>
    static
    typename std::enable_if<NodalDOF != 0,
        decltype(std::declval<Mat>().template block<NIntPts, NIntPts>(0u, 0u))
    >::type
    blockShpShp(Mat& mat,
                unsigned top, unsigned left, unsigned width, unsigned height)
    {
        assert(width==NIntPts && height==NIntPts);
        (void) width; (void) height;
        return mat.template block<NIntPts, NIntPts>(top, left);
    }
    // block gauss x gauss, dynamic-size non-const matrix
    template<typename Mat>
    static
    typename std::enable_if<NodalDOF == 0,
        decltype(std::declval<Mat>().block(0u, 0u, 0u, 0u))
    >::type
    blockShpShp(Mat& mat,
                unsigned top, unsigned left, unsigned width, unsigned height)
    {
        assert(width == height);
        return mat.block(top, left, width, height);
    }

    // block gauss x 1, fixed-size const vector
    template<typename Vec>
    static
    typename std::enable_if<NodalDOF != 0,
        decltype(std::declval<const Vec>().template block<NIntPts, 1>(0u, 0u))
    >::type
    blockShp(Vec const& vec, unsigned top, unsigned height)
    {
        assert(height==NIntPts);
        (void) height;
        return vec.template block<NIntPts, 1>(top, 0);
    }
    // block gauss x 1, dynamic-size const vector
    template<typename Vec>
    static
    typename std::enable_if<NodalDOF == 0,
        decltype(std::declval<const Vec>().block(0u, 0u, 0u, 0u))
    >::type
    blockShp(Vec const& vec, unsigned top, unsigned height)
    {
        return vec.block(top, 0, height, 1);
    }
    // block gauss x 1, fixed-size non-const vector
    template<typename Vec>
    static
    typename std::enable_if<NodalDOF != 0,
        decltype(std::declval<Vec>().template block<NIntPts, 1>(0u, 0u))
    >::type
    blockShp(Vec& vec, unsigned top, unsigned height)
    {
        assert(height==NIntPts);
        (void) height;
        return vec.template block<NIntPts, 1>(top, 0);
    }
    // block gauss x 1, dynamic-size non-const vector
    template<typename Vec>
    static
    typename std::enable_if<NodalDOF == 0,
        decltype(std::declval<Vec>().block(0u, 0u, 0u, 0u))
    >::type
    blockShp(Vec& vec, unsigned top, unsigned height)
    {
        return vec.block(top, 0, height, 1);
    }

private:
    static const unsigned _Dim = Dim;
    static const unsigned _NIntPts = NIntPts;
};

#ifndef EIGEN_DYNAMIC_SHAPE_MATRICES

template<typename ShpPol, unsigned NIntPts, unsigned NodalDOF, unsigned Dim>
using DataTraits = DataTraitsFixed<ShpPol, NIntPts, NodalDOF, Dim>;

static_assert(EIGEN_DYNAMIC_SHAPE_MATRICES_FLAG == 0, "inconsistent use of macros");
#else

template<typename ShpPol, unsigned NIntPts, unsigned NodalDOF, unsigned Dim>
using DataTraits = DataTraitsFixed<ShapeMatrixPolicyType<void, 0>, 0, 0, 0>;

static_assert(EIGEN_DYNAMIC_SHAPE_MATRICES_FLAG == 1, "inconsistent use of macros");
#endif


template<typename Traits>
class TESFEMReactionAdaptor;

template<typename Traits>
class TESFEMReactionAdaptorAdsorption;

template<typename Traits>
class TESFEMReactionAdaptorInert;

template<typename Traits>
class TESFEMReactionAdaptorSinusoidal;

template<typename Traits>
class TESFEMReactionAdaptorCaOH2;


template<typename Traits>
class LADataNoTpl
{
public:
    void assembleIntegrationPoint(
            unsigned integration_point,
            std::vector<double> const& localX,
            typename Traits::ShapeMatrices::ShapeType const& smN,
            typename Traits::ShapeMatrices::DxShapeType const& smDNdx,
            typename Traits::ShapeMatrices::JacobianType const& smJ,
            const double smDetJ,
            const double weight
            );

    void init(const unsigned num_int_pts, const unsigned dimension);

    void preEachAssemble();
    void postEachAssemble(typename Traits::LocalMatrix& local_M,
                          typename Traits::LocalMatrix& local_K,
                          typename Traits::LocalVector& local_b);

    std::vector<double> const&
    getIntegrationPointValues(SecondaryVariables var, std::vector<double>& cache) const;

private:
    Eigen::Matrix3d getMassCoeffMatrix(const unsigned int_pt);
    typename Traits::LaplaceMatrix getLaplaceCoeffMatrix(const unsigned int_pt, const unsigned dim);
    Eigen::Matrix3d getAdvectionCoeffMatrix(const unsigned int_pt);
    Eigen::Matrix3d getContentCoeffMatrix(const unsigned int_pt);
    Eigen::Vector3d getRHSCoeffVector(const unsigned int_pt);

    void preEachAssembleIntegrationPoint(
            const unsigned int_pt,
            std::vector<double> const& localX,
            typename Traits::ShapeMatrices::ShapeType const& smN,
            typename Traits::ShapeMatrices::DxShapeType const& smDNdx,
            typename Traits::ShapeMatrices::JacobianType const& smJ,
            const double smDetJ
            );

    void initReaction(const unsigned int_pt);

    // TODO data members except local matrices are independent of any template parameter
    // they can be moved to a separate non-template struct for better decoupling of the
    // reaction adaptor.
    // Maybe the reaction adaptor does not even need direct access to those members!

    AssemblyParams const* _AP;

    // integration point quantities
    std::vector<double> _solid_density;
    std::vector<double> _reaction_rate; // dC/dt * _rho_SR_dry
    std::vector<std::vector<double> > _velocity; // vector of velocities for each integration point


    // integration point values of unknowns -- temporary storage
    double _p = std::numeric_limits<double>::quiet_NaN(); // gas pressure
    double _T = std::numeric_limits<double>::quiet_NaN(); // temperature
    double _vapour_mass_fraction = std::numeric_limits<double>::quiet_NaN(); // fluid mass fraction of the second component

    // temporary storage for some properties
    // values do not change during the assembly of one integration point
    double _rho_GR = std::numeric_limits<double>::quiet_NaN();
    double _p_V    = std::numeric_limits<double>::quiet_NaN(); // vapour partial pressure
    double _qR     = std::numeric_limits<double>::quiet_NaN();  // reaction rate, use this in assembly!!!

    // TODO: entirely omit local matrices
    typename Traits::LocalMatrix _Mas;
    typename Traits::LocalMatrix _Lap_Adv_Cnt;
    typename Traits::LocalVector _rhs;

    std::unique_ptr<TESFEMReactionAdaptor<Traits> > _reaction_adaptor;

    // variables at previous timestep
    std::vector<double> _solid_density_prev_ts;
    std::vector<double> _reaction_rate_prev_ts; // could also be calculated from _solid_density_prev_ts


    template <typename, typename, typename, typename, unsigned>
    friend class LocalAssemblerData;

    friend class TESFEMReactionAdaptor<Traits>;
    friend class TESFEMReactionAdaptorAdsorption<Traits>;
    friend class TESFEMReactionAdaptorInert<Traits>;
    friend class TESFEMReactionAdaptorSinusoidal<Traits>;
    friend class TESFEMReactionAdaptorCaOH2<Traits>;
};


template <typename Vec>
void
ogs5OutVec(const Vec& vec);

template <typename Mat>
void
ogs5OutMat(const Mat& mat);


} // namespace TES

} // namespace ProcessLib


// tricking cmake dependency checker
#include "TESFEM-data-impl-incl.h"

#endif // PROCESS_LIB_TES_FEM_NOTPL_H_
