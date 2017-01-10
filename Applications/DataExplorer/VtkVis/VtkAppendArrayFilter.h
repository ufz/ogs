/**
 * \file
 * \author Karsten Rink
 * \date   2010-02-09
 * \brief  Definition of the VtkAppendArrayFilter class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

// ** INCLUDES **
#include "VtkAlgorithmProperties.h"

#include <vtkUnstructuredGridAlgorithm.h>

#include <vector>

/**
 * \brief
 */
class VtkAppendArrayFilter : public vtkUnstructuredGridAlgorithm, public VtkAlgorithmProperties
{
public:
    /// @brief Create new objects with New() because of VTKs object reference counting.
    static VtkAppendArrayFilter* New();

    vtkTypeMacro(VtkAppendArrayFilter, vtkUnstructuredGridAlgorithm);

    /// @brief Prints the mesh data to an output stream.
    void PrintSelf(ostream& os, vtkIndent indent) override;

    /// @brief Sets user properties.
    void SetUserProperty(QString name, QVariant value) override
    {
        Q_UNUSED(name);
        Q_UNUSED(value);
    }

    void SetArray(const std::string &array_name, const std::vector<double> &selection);

protected:
    VtkAppendArrayFilter();
    ~VtkAppendArrayFilter();

    /// @brief The filter logic.
    int RequestData(vtkInformation* request,
                    vtkInformationVector** inputVector,
                    vtkInformationVector* outputVector) override;

private:
    std::vector<double> _array;
    std::string _array_name;
};
