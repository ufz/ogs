#ifndef PROCESS_LIB_DensityDrivenFlowMATERIALMODEL_H_
#define PROCESS_LIB_DensityDrivenFlowMATERIALMODEL_H_

#include <math.h>

namespace ProcessLib
{
namespace DensityDrivenFlow
{

static double DensityWater_T(double density0, double temperature, double temperature0, double beta)
{
   return density0*(1 - beta*(temperature - temperature0));
}

static double Viscosity(double viscosity0, double temperature, double temperature_con, double temperature0)
{
    return viscosity0*exp(-(temperature - temperature0)/temperature_con); 
}

}  // DensityDrivenFlow

}  // ProcessLib


#endif  // PROCESS_LIB_DensityDrivenFlowMATERIALMODEL_H_
