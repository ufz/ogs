/*!
   \file  LinearSolidDensityModel.h
   \brief Declaration of class LinearSolidDensityModel for solid density
          depending on a variable linearly.
   
   \author Wenqing Wang
   \date Jan 2015

   \copyright
    Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
               Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license
*/
#ifndef LINEAR_SOLID_DENSITY_MODEL_H_
#define LINEAR_SOLID_DENSITY_MODEL_H_

#include <string>

#include "MaterialLib/DensityType.h"

namespace MaterialLib
{

/// Solid model depends on a variable
class LinearSolidDensityModel
{
    public:
        /*!
             \param var1      Reference variable 1.
             \param density1  Density at variable 1.
             \param var2  Reference variable 2.
             \param density2  Density at variable 2.
             
        */
        LinearSolidDensityModel(const double var1, const double density1,
                                const double var2, const double density2)
                     :_var1(var1), _density1(density1),
                      _tangential( (density2-density1)/(var2-var1) ) 
        {
        }
        
        /// Get desity model name.
        std::string getName() const
        {
           return "Linear density of solid" ;
        }

        DensityType getType() const
        {
            return DensityType::SOLID_LINEAR;
        }
        /// Get density value
        /// \param var Variable
        double getValue(const double var) const
        {
           return _tangential * (var - _var1) + _density1;  
        }
    private:
       double _var1;       ///<  Reference variable 1.
       double _density1;   ///<  Reference density 1.
       double _tangential; ///< Tangential of the linear curve.
};

} // end namespace
#endif

