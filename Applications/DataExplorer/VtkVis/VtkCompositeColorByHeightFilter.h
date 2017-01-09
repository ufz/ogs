/**
 * \file
 * \author Karsten Rink
 * \date   2010-11-01
 * \brief  Definition of the VtkCompositeColorByHeightFilter class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef VTKCOMPOSITECOLORBYHEIGHTFILTER_H
#define VTKCOMPOSITECOLORBYHEIGHTFILTER_H

#include "VtkCompositeFilter.h"

class vtkSphereSource;

/// @brief This filter colors the input by the points z-value.
class VtkCompositeColorByHeightFilter : public VtkCompositeFilter
{
public:
    VtkCompositeColorByHeightFilter(vtkAlgorithm* inputAlgorithm);
    virtual ~VtkCompositeColorByHeightFilter() {}

    virtual void init();

    virtual void SetUserProperty(QString name, QVariant value);

protected:
};

#endif // VTKCOMPOSITECOLORBYHEIGHTFILTER_H
