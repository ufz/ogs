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

#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "ProcessLib/LocalAssemblerTraits.h"
#include "ProcessLib/VariableTransformation.h"

#include "TESAssemblyParams.h"

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


class TESFEMReactionAdaptor;

struct TESLocalAssemblerData
{
    AssemblyParams const* ap;

    // integration point quantities
    std::vector<double> solid_density;
    std::vector<double> reaction_rate; // dC/dt * rho_SR_dry
    std::vector<std::vector<double> > velocity; // vector of velocities for each integration point


    // integration point values of unknowns -- temporary storage
    double p = std::numeric_limits<double>::quiet_NaN(); // gas pressure
    double T = std::numeric_limits<double>::quiet_NaN(); // temperature
    double vapour_mass_fraction = std::numeric_limits<double>::quiet_NaN(); // fluid mass fraction of the second component

    // temporary storage for some properties
    // values do not change during the assembly of one integration point
    double rho_GR = std::numeric_limits<double>::quiet_NaN();
    double p_V    = std::numeric_limits<double>::quiet_NaN(); // vapour partial pressure
    double qR     = std::numeric_limits<double>::quiet_NaN();  // reaction rate, use this in assembly!!!

    std::unique_ptr<TESFEMReactionAdaptor> reaction_adaptor;

    // variables at previous timestep
    std::vector<double> solid_density_prev_ts;
    std::vector<double> reaction_rate_prev_ts; // could also be calculated from solid_density_prev_ts
};


template<typename Traits>
class TESLocalAssemblerInner
{
public:
    void assembleIntegrationPoint(
            unsigned integration_point,
            std::vector<double> const& localX,
            typename Traits::ShapeMatrices::ShapeType const& smN,
            typename Traits::ShapeMatrices::DxShapeType const& smDNdx,
            typename Traits::ShapeMatrices::JacobianType const& smJ,
            const double smDetJ,
            const double weight,
            typename Traits::LocalMatrix& local_M,
            typename Traits::LocalMatrix& local_K,
            typename Traits::LocalVector& local_b
            );

    void init(const unsigned num_int_pts, const unsigned dimension);

    void preEachAssemble();

    std::vector<double> const&
    getIntegrationPointValues(SecondaryVariables var, std::vector<double>& cache) const;

    // TODO pass to constructor
    void setAssemblyParameters(AssemblyParams const& ap) { _d.ap = &ap; }
    // TODO better encapsulation
    AssemblyParams const& getAssemblyParameters() const { return * _d.ap; }
    TESFEMReactionAdaptor const& getReactionAdaptor() const {
        return *_d.reaction_adaptor;
    }
    TESFEMReactionAdaptor& getReactionAdaptor() {
        return *_d.reaction_adaptor;
    }

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

    TESLocalAssemblerData _d;
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
#include "TESLocalAssemblerInner-impl-incl.h"

#endif // PROCESS_LIB_TES_FEM_NOTPL_H_
