/**
 * \file
 * \author Karsten Rink
 * \date   2011-09-29
 * \brief  Definition of the VtkVisPointSetItem class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef VTKVISPOINTSETITEM_H
#define VTKVISPOINTSETITEM_H

// ** INCLUDES **
#include "VtkVisPipelineItem.h"

class vtkAlgorithm;
class vtkDataSetAttributes;
class vtkPointSet;
class vtkProp3D;
class vtkRenderer;
class vtkTransformFilter;
class QVtkDataSetMapper;

class VtkAlgorithmProperties;
class VtkCompositeFilter;

/**
 * \brief An item in the VtkVisPipeline containing a point set object to be visualized.
 *
 * Any VTK point set object (i.e. vtkUnstructuredGrid- and vtkPolyDataAlgorithm-objects)
 * are represented by a VtkVisPointSetItem to be assigned a mapper, an actor and its
 * visualization properties (colour, scalar values, etc.).
 * \sa VtkVisPipelineItem
 */
class VtkVisPointSetItem : public VtkVisPipelineItem
{

public:
    /// @brief Constructor for a source/filter object.
    VtkVisPointSetItem(vtkAlgorithm* algorithm,
                       TreeItem* parentItem,
                       const QList<QVariant> data = QList<QVariant>());

    /// @brief Constructor for composite filter
    VtkVisPointSetItem(VtkCompositeFilter* compositeFilter, TreeItem* parentItem,
                       const QList<QVariant> data = QList<QVariant>());

    ~VtkVisPointSetItem();

    /// @brief Gets the last selected attribute.
    const QString GetActiveAttribute() const;

    /// @brief Get the scalar range for the active attribute
    void GetRangeForActiveAttribute(double range[2]) const;

    /// @brief Initializes vtkMapper and vtkActor necessary for visualization of
    /// the item and sets the item's properties.
    void Initialize(vtkRenderer* renderer);

    vtkAlgorithm* transformFilter() const;

    /// @brief Sets the selected attribute array for the visualisation of the data set.
    void SetActiveAttribute(const QString& name);

    /// @brief Scales the data in visualisation-space.
    void setScale(double x, double y, double z) const;

    /// @brief Translates the item in visualisation-space.
    void setTranslation(double x, double y, double z) const;

    /// @brief Enables / disables backface culling.
    void setBackfaceCulling(bool enable) const;

protected:
    QVtkDataSetMapper* _mapper;
    vtkTransformFilter* _transformFilter;
    bool _onPointData;
    std::string _activeArrayName;

    /// Selects the appropriate VTK-Writer object and writes the object to a file with the given name.
    virtual int callVTKWriter(vtkAlgorithm* algorithm, const std::string &filename) const;

    void SetScalarVisibility(bool on);

    /// @brief Sets pre-set properties on vtkActor and on vtkMapper
    void setVtkProperties(VtkAlgorithmProperties* vtkProps);

private:
    /// Checks if the selected attribute actually exists for the data set
    bool activeAttributeExists(vtkDataSetAttributes* data, std::string& name);

};

#endif // VTKVISPOINTSETITEM_H
