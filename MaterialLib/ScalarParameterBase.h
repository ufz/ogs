/*!
   \file  ScalarParameterBase.h
   \brief Declaration of root class for scalar material parameter.

   \author Wenqing Wang
   \date Jan 2015

   \copyright
    Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
               Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license
*/
#ifndef SCALAR_PARAMETER_BASE_H_
#define SCALAR_PARAMETER_BASE_H_

namespace MaterialLib
{
/// Base scalar parameter class.
template<typename T_PARAMETER_TYPE> class ScalarParameterBase
{
    public:
        void setType(const T_PARAMETER_TYPE type )
        {
            _type = type;
        }

        T_PARAMETER_TYPE getType() const
        {
            return _type;
        }

        double getValue()
        {
            return 0.;
        }

    private:
        T_PARAMETER_TYPE _type = static_cast<T_PARAMETER_TYPE>(0);
};

} // end namespace
#endif
