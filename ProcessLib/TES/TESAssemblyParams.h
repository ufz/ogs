#ifndef PROCESS_LIB_TESPROCESS_NOTPL_H_
#define PROCESS_LIB_TESPROCESS_NOTPL_H_


#include <Eigen/Sparse>
#include <Eigen/Eigen>

#include "MaterialsLib/adsorption/reaction.h"

#include "ProcessLib/VariableTransformation.h"


namespace ProcessLib
{

namespace TES
{

const unsigned NODAL_DOF = 3;

const double M_N2  = 0.028013;
const double M_H2O = 0.018016;

struct AssemblyParams
{
    Trafo trafo_p;
    Trafo trafo_T;
    Trafo trafo_x;

    std::unique_ptr<Ads::Reaction> react_sys;

    double fluid_specific_heat_source = std::numeric_limits<double>::quiet_NaN();
    double cpG = std::numeric_limits<double>::quiet_NaN(); // specific isobaric fluid heat capacity

    Eigen::MatrixXd solid_perm_tensor =
            Eigen::MatrixXd::Constant(3, 3, std::numeric_limits<double>::quiet_NaN()); // TODO get dimensions
    double solid_specific_heat_source = std::numeric_limits<double>::quiet_NaN();
    double solid_heat_cond = std::numeric_limits<double>::quiet_NaN();
    double cpS = std::numeric_limits<double>::quiet_NaN();    // specific isobaric solid heat capacity

    double tortuosity = std::numeric_limits<double>::quiet_NaN();
    double diffusion_coefficient_component = std::numeric_limits<double>::quiet_NaN(); // ???

    double poro = std::numeric_limits<double>::quiet_NaN();

    double rho_SR_dry = std::numeric_limits<double>::quiet_NaN();

    const double M_inert = M_N2; // N2
    const double M_react = M_H2O;

    // TODO unify variable names
    double initial_solid_density = std::numeric_limits<double>::quiet_NaN();

    double       delta_t = std::numeric_limits<double>::quiet_NaN();
    unsigned     iteration_in_current_timestep = 0;

    bool output_element_matrices = false;

    unsigned number_of_try_of_iteration = 0;
    double   current_time = std::numeric_limits<double>::quiet_NaN();
};


class TESProcessInterface
{
public:
    AssemblyParams const& getAssemblyParams() const {
        return _assembly_params;
    }

    virtual ~TESProcessInterface() = default;

protected:
    AssemblyParams _assembly_params;
};

} // namespace TES

} // namespace ProcessLib

#endif // PROCESS_LIB_TESPROCESS_NOTPL_H_
