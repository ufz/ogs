/*!
   \file  ParameterBase.h
   \brief Declaration of root class for material parameter classes,
          which holds parameter type only,

   \author Wenqing Wang
   \date Jan 2015

   \copyright
    Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
               Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license
*/
#ifndef PARAMETER_BASE_H_
#define PARAMETER_BASE_H_

namespace MaterialLib
{
/// Root class for material parameter classes.
template<typename T_PARAMETER_TYPE> class ParameterBase
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

    private:
        T_PARAMETER_TYPE _type = static_cast<T_PARAMETER_TYPE>(0);
};

} // end namespace
#endif
