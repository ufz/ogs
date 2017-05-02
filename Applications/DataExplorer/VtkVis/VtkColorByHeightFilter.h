/**
 * \file
 * \author Karsten Rink
 * \date   2010-04-21
 * \brief  Definition of the VtkColorByHeightFilter class.
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

#include <vtkPolyDataAlgorithm.h>
#include <vtkType.h>

class VtkColorLookupTable;

/**
 * \brief VTK filter object for colouring vtkPolyData objects based on z-coordinates.
 *
 * This filter class is basically a container for a ColorLookupTable. In fact, you can get the underlying
 * ColorLookupTable using the method GetColorLookupTable(). Using this method allows the user to set a number
 * of properties on that lookup table such as interpolation method, the range of values over which the lookup
 * table is calculated and so on.
 * If no range boundaries are explicitly set, the minimum and maximum height value will be calculated from
 * the data and set as minimum and maximum value for the lookup table.
 * ColorLookupTable must be deleted manually.
 * \see VtkCompositeColorByHeightFilter::init() for example usage.
 */
class VtkColorByHeightFilter : public vtkPolyDataAlgorithm, public VtkAlgorithmProperties
{
public:
    /// @brief Create new objects with New() because of VTKs object reference counting.
    static VtkColorByHeightFilter* New();

    vtkTypeMacro(VtkColorByHeightFilter, vtkPolyDataAlgorithm);

    /// @brief Prints the mesh data to an output stream.
    void PrintSelf(ostream& os, vtkIndent indent) override;

    /// @brief Returns the underlying colour look up table object.
    vtkGetObjectMacro(ColorLookupTable,VtkColorLookupTable);

    /// @brief This filter gets updated when the color look-up table was modified.
    vtkMTimeType GetMTime() override;

    /// @brief Sets user properties.
    void SetUserProperty(QString name, QVariant value) override
    {
        Q_UNUSED(name);
        Q_UNUSED(value);
    }

    /// @brief Sets the boundaries for the color look-up table.
    void SetTableRange(double min, double max);

    /// @brief Sets the scaling of the color look-up table boundaries.
    /// This is used in VtkVisTabWidget when a parent filter is scaled.
    void SetTableRangeScaling(double scale);

protected:
    VtkColorByHeightFilter();
    ~VtkColorByHeightFilter() override;

    /// @brief The filter logic.
    int RequestData(vtkInformation* request,
                    vtkInformationVector** inputVector,
                    vtkInformationVector* outputVector) override;

    /// @brief Calculates the color lookup table based on set parameters.
    VtkColorLookupTable* BuildColorTable();

    VtkColorLookupTable* ColorLookupTable;

    double _tableRange[2];
    double _tableRangeScaling;
};
