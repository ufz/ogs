/**
 * \file
 * \author Lars Bilke
 * \date   2010-10-20
 * \brief  Definition of the VtkFilterFactory class.
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <QString>
#include <QVector>
#include <vtkSetGet.h>

class VtkCompositeFilter;
class vtkAlgorithm;
struct VtkFilterInfo;

/// @brief Creates registered filter objects by name.
class VtkFilterFactory
{
public:
    //VtkFilterFactory();
    //virtual ~VtkFilterFactory();

    /// @brief Returns all registered filters.
    /// New VtkCompositeFilter or filter inherited from VtkAlgorithmProperties
    /// must be registered here.
    static const QVector<VtkFilterInfo> GetFilterList();

    /// @brief Creates a composite filter by name.
    static VtkCompositeFilter* CreateCompositeFilter(QString type, vtkAlgorithm* inputAlgorithm);

    /// @brief Creates a normal filter name.
    static vtkAlgorithm* CreateSimpleFilter(QString type);
};

/// @brief Holds meta information about a filter
struct VtkFilterInfo
{
    /// @brief Constructor.
    /// @param name The name of the filter (the class name)
    /// @param readableName
    /// @param description A short description of what the filter does
    /// @param inputDataObjectType The input data type (see OutputDataObjectTypeAsString())
    /// @param outputDataObjectType The output data type (see OutputDataObjectTypeAsString())
    VtkFilterInfo(QString name, QString readableName, QString description,
                  int inputDataObjectType, int outputDataObjectType)
    {
        this->name = name;
        this->readableName = readableName;
        this->description = description;
        this->inputDataObjectType = inputDataObjectType;
        this->outputDataObjectType = outputDataObjectType;
    }

    /// @brief Default constructor.
    VtkFilterInfo()
    {
        this->name = QString();
        this->readableName = QString();
        this->description = QString();
        this->inputDataObjectType = -1;
        this->outputDataObjectType = -1;
    }

    /// @brief Returns the data type as a string.
    QString OutputDataObjectTypeAsString() const
    {
        switch (outputDataObjectType)
        {
        case VTK_POLY_DATA: return QString("vtkPolyData");
        case VTK_STRUCTURED_POINTS: return QString("vtkStructuredPoints");
        case VTK_STRUCTURED_GRID: return QString("vtkStructuredGrid");
        case VTK_RECTILINEAR_GRID: return QString("vtkRectilinearGrid");
        case VTK_UNSTRUCTURED_GRID: return QString("vtkUnstructuredGrid");
        case VTK_IMAGE_DATA: return QString("vtkImageData");
        case VTK_DATA_SET: return QString("vtkDataSet");
        default: return QString("Data type not defined!");
        }
    }

    QString name;
    QString readableName;
    QString description;
    int inputDataObjectType;
    int outputDataObjectType;
};
