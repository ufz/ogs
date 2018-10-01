/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

/**
 * \file BHEAbstract.h
 * 2014/06/04 HS inital implementation
 * borehole heat exchanger abstract class
 *
 * 1) Diersch_2011_CG
 * Two very important references to understand this class implementations are:
 * H.-J.G. Diersch, D. Bauer, W. Heidemann, W. R黨aak, P. Sch鋞zl,
 * Finite element modeling of borehole heat exchanger systems:
 * Part 1. Fundamentals, Computers & Geosciences,
 * Volume 37, Issue 8, August 2011, Pages 1122-1135, ISSN 0098-3004,
 * http://dx.doi.org/10.1016/j.cageo.2010.08.003.
 *
 * 2) FEFLOW_2014_Springer
 * FEFLOW: Finite Element Modeling of Flow, Mass and Heat Transport in Porous
 * and Fractured Media Diersch, Hans-Joerg, 2014, XXXV, 996 p, Springer.
 *
 */

#pragma once

#include <Eigen/Eigen>
#include <map>
#include "GeoLib/Polyline.h"
#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"
#include "ProcessLib/Utils/ProcessUtils.h"
#include "boost/math/constants/constants.hpp"

namespace ProcessLib
{
namespace HeatTransportBHE
{
namespace BHE  // namespace of borehole heat exchanger
{
/**
 * list the types of borehole heat exchanger
 */
enum class BHE_TYPE
{
    TYPE_2U,   // two u-tube borehole heat exchanger
    TYPE_1U,   // one u-tube borehole heat exchanger
    TYPE_CXC,  // coaxial pipe with annualar inlet
    TYPE_CXA   // coaxial pipe with centreed inlet
};

enum class BHE_BOUNDARY_TYPE
{
    FIXED_INFLOW_TEMP_BOUNDARY,
    FIXED_INFLOW_TEMP_CURVE_BOUNDARY,
    POWER_IN_WATT_BOUNDARY,
    POWER_IN_WATT_CURVE_FIXED_DT_BOUNDARY,
    BUILDING_POWER_IN_WATT_CURVE_FIXED_DT_BOUNDARY,
    BUILDING_POWER_IN_WATT_CURVE_FIXED_FLOW_RATE_BOUNDARY,
    POWER_IN_WATT_CURVE_FIXED_FLOW_RATE_BOUNDARY,
    FIXED_TEMP_DIFF_BOUNDARY
};

/**
 * discharge type of the 2U BHE
 */
enum class BHE_DISCHARGE_TYPE
{
    BHE_DISCHARGE_TYPE_PARALLEL,  // parallel discharge
    BHE_DISCHARGE_TYPE_SERIAL     // serial discharge
};

class BHEAbstract
{
public:
    struct BoreholeGeometry
    {
        /**
         * length/depth of the BHE
         * unit is m
         */
        double length;

        /**
         * diameter of the BHE
         * unit is m
         */
        double diameter;
    };

    struct PipeParameters
    {
        /**
         * radius of the pipline inner side
         * unit is m
         */
        double r_inner;

        /**
         * radius of the pipline outer side
         * unit is m
         */
        double r_outer;

        /**
         * pipe-in wall thickness
         * unit is m
         */
        double b_in;

        /**
         * pipe-out wall thickness
         * unit is m
         */
        double b_out;

        /**
         * thermal conductivity of the pipe wall
         * unit is kg m sec^-3 K^-1
         */

        double lambda_p;

        /**
         * thermal conductivity of the inner pipe wall
         * unit is kg m sec^-3 K^-1
         */
        double lambda_p_i;

        /**
         * thermal conductivity of the outer pipe wall
         * unit is kg m sec^-3 K^-1
         */
        double lambda_p_o;
    };

    struct RefrigerantParameters
    {
        /**
         * dynamics viscosity of the refrigerant
         * unit is kg m-1 sec-1
         */
        double mu_r;

        /**
         * density of the refrigerant
         * unit is kg m-3
         */
        double rho_r;

        /**
         * thermal conductivity of the refrigerant
         * unit is kg m sec^-3 K^-1
         */
        double lambda_r;

        /**
         * specific heat capacity of the refrigerant
         * unit is m^2 sec^-2 K^-1
         */
        double heat_cap_r;

        /**
         * longitudinal dispersivity of the
         * referigerant flow in the pipeline
         */
        double alpha_L;
    };

    struct GroutParameters
    {
        /**
         * density of the grout
         * unit is kg m-3
         */
        double rho_g;

        /**
         * porosity of the grout
         * unit is [-]
         */
        double porosity_g;

        /**
         * specific heat capacity of the grout
         * unit is m^2 sec^-2 K^-1
         */
        double heat_cap_g;

        /**
         * thermal conductivity of the grout
         * unit is kg m sec^-3 K^-1
         */
        double lambda_g;
    };

    struct ExternallyDefinedRaRb
    {
        /**
         * whether or not using external given borehole thermal resistance
         * values Ra, Rb
         */
        bool use_extern_Ra_Rb;

        /**
         * external given borehole internal thermal resistance value
         */
        double ext_Ra;

        /**
         * external given borehole thermal resistance value
         */
        double ext_Rb;
    };

    struct ExternallyDefinedThermalResistances
    {
        /**
         * whether or not using user defined borehole thermal resistance Rfig,
         * Rfog, Rgg, Rgs
         */
        bool if_use_defined_therm_resis;

        /**
         * external given borehole internal thermal resistance value
         */
        double ext_Rfig;

        /**
         * external given borehole internal thermal resistance value
         */
        double ext_Rfog;

        /**
         * external given borehole internal thermal resistance value
         */
        double ext_Rgg1;

        /**
         * external given borehole internal thermal resistance value
         */
        double ext_Rgg2;

        /**
         * external given borehole internal thermal resistance value
         */
        double ext_Rgs;
    };

    /**
     * constructor
     */
    BHEAbstract(
        const std::string name_,
        const BHE_TYPE bhe_type_,
        BoreholeGeometry borehole_geometry_,
        PipeParameters pipe_param_,
        RefrigerantParameters refrigerant_param_,
        GroutParameters grout_param_,
        ExternallyDefinedRaRb extern_Ra_Rb_,
        ExternallyDefinedThermalResistances extern_def_thermal_resistances_,
        std::map<std::string,
                 std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
            bhe_curves_,
        BHE_BOUNDARY_TYPE my_bound_type =
            BHE_BOUNDARY_TYPE::FIXED_INFLOW_TEMP_BOUNDARY,
        bool if_use_ext_Ra_Rb = false,
        bool user_defined_R_vals = false,
        bool if_flowrate_curve = false)
        : name(name_),
          bhe_type(bhe_type_),
          boundary_type(my_bound_type),
          borehole_geometry(borehole_geometry_),
          pipe_param(pipe_param_),
          refrigerant_param(refrigerant_param_),
          grout_param(grout_param_),
          extern_Ra_Rb(extern_Ra_Rb_),
          extern_def_thermal_resistances(extern_def_thermal_resistances_),
          bhe_curves(bhe_curves_),
          use_flowrate_curve(if_flowrate_curve),
          if_use_ext_Ra_Rb(if_use_ext_Ra_Rb),
          user_defined_R_vals(user_defined_R_vals),
          PI(boost::math::constants::pi<double>()){};

    /**
     * destructor
     */
    virtual ~BHEAbstract() = default;

    /**
     * return the number of unknowns needed for this BHE
     * abstract function, need to be realized.
     */
    virtual std::size_t getNumUnknowns() const = 0;

    /**
     * initialization calcultion,
     * need to be overwritten.
     */
    virtual void initialize()
    {
        calcPipeFlowVelocity();
        calcRenoldsNum();
        calcPrandtlNum();
        calcNusseltNum();
        calcThermalResistances();
        calcHeatTransferCoefficients();
    };

    /**
     * update all parameters based on the new flow rate
     * not necessarily needs to be overwritten.
     */
    virtual void updateFlowRate(double new_flow_rate)
    {
        Q_r = new_flow_rate;
        calcPipeFlowVelocity();
        calcRenoldsNum();
        calcPrandtlNum();
        calcNusseltNum();
        calcThermalResistances();
        calcHeatTransferCoefficients();
    };

    virtual void updateFlowRateFromCurve(double current_time) = 0;

    /**
     * thermal resistance calculation,
     * need to be overwritten.
     */
    virtual void calcThermalResistances() = 0;

    /**
     * Nusselt number calculation,
     * need to be overwritten.
     */
    virtual void calcNusseltNum() = 0;

    /**
     * Renolds number calculation,
     * need to be overwritten.
     */
    virtual void calcRenoldsNum() = 0;

    /**
     * Prandtl number calculation,
     * need to be overwritten.
     */
    virtual void calcPrandtlNum() = 0;

    /**
     * flow velocity inside the pipeline
     * need to be overwritten.
     */
    virtual void calcPipeFlowVelocity() = 0;

    /**
     * heat transfer coefficient,
     * need to be overwritten.
     */
    virtual void calcHeatTransferCoefficients() = 0;

    /**
     * return the coeff of mass matrix,
     * depending on the index of unknown.
     */
    virtual double getMassCoeff(std::size_t idx_unknown) const = 0;

    /**
     * return the coeff of laplace matrix,
     * depending on the index of unknown.
     */
    virtual void getLaplaceMatrix(std::size_t idx_unknown,
                                  Eigen::MatrixXd& mat_laplace) const = 0;

    /**
     * return the coeff of advection matrix,
     * depending on the index of unknown.
     */
    virtual void getAdvectionVector(std::size_t idx_unknown,
                                    Eigen::VectorXd& vec_advection) const = 0;

    /**
     * return the coeff of boundary heat exchange matrix,
     * depending on the index of unknown.
     */
    virtual double getBoundaryHeatExchangeCoeff(
        std::size_t idx_unknown) const = 0;

    /**
     * return the inflow temperature based on outflow temperature and fixed
     * power.
     */
    virtual double getTinByTout(double T_in, double current_time) = 0;

    /**
     * return the _R_matrix, _R_pi_s_matrix, _R_s_matrix
     */
    virtual void setRMatrices(
        const int idx_bhe_unknowns, const int NumNodes,
        Eigen::MatrixXd& matBHE_loc_R,
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>&
            R_matrix,
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>&
            R_pi_s_matrix,
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>&
            R_s_matrix) const = 0;

public:
    /**
     * name of the borehole heat exchanger
     */
    const std::string name;

    /**
     * the type of this BHE
     */
    const BHE_TYPE bhe_type;

    /**
     * the type of the boundary condition on this BHE
     */
    const BHE_BOUNDARY_TYPE boundary_type;

    /**
     * geometry of the borehole
     */
    BoreholeGeometry const borehole_geometry;

    /**
     * geometry of the pipes in the borehole
     */
    PipeParameters const pipe_param;

    /**
     * parameters of the refrigerant
     */
    RefrigerantParameters const refrigerant_param;

    /**
     * parameters of the grout
     */
    GroutParameters const grout_param;

    /**
     * Ra Rb values defined by the user externally
     */
    ExternallyDefinedRaRb const extern_Ra_Rb;

    /**
     * thermal resistance values defined by the user externally
     */
    ExternallyDefinedThermalResistances const extern_def_thermal_resistances;

    /**
     * power extracted from or injected into the BHE
     * unit is Watt
     * if value positive, then injecting power
     * if value negative, then extracting power
     */
    double power_in_watt_val;

    /**
     * temperature difference between inflow and
     * outflow pipelines
     */
    double delta_T_val;

    /**
     * threshold Q value for switching off the BHE
     * when using the Q_curve_fixed_dT B.C.
     */
    double threshold;

    /**
     * map strucutre that contains all the curves related to this BHE
     */
    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
        bhe_curves;

    /**
     * power in watt curve
     */
    MathLib::PiecewiseLinearInterpolation* power_in_watt_curve;

    /**
     * heating COP curve
     */
    MathLib::PiecewiseLinearInterpolation* heating_cop_curve;

    /**
     * cooling COP curve
     */
    MathLib::PiecewiseLinearInterpolation* cooling_cop_curve;

    /**
     * refrigerant flow rate curve
     */
    MathLib::PiecewiseLinearInterpolation* flowrate_curve;

    /**
     * inflow temperature curve
     */
    MathLib::PiecewiseLinearInterpolation* inflow_temperature_curve;

    /**
     * use refrigerant flow rate curve
     */
    bool use_flowrate_curve;

    /**
     * whether external Ra an Rb values are supplied by the user
     */
    const bool if_use_ext_Ra_Rb;

    /**
     * whether external R values are supplied by the user
     */
    const bool user_defined_R_vals;

protected:
    /**
     * total refrigerant flow discharge of BHE
     * unit is m^3/sec
     */
    double Q_r;

    const double PI;
};
}  // end of namespace BHE
}  // end of namespace HeatTransportBHE
}  // end of namespace ProcessLib
