/**
 * \file
 * \author Karsten Rink
 * \date   2011-09-29
 * \brief  Definition of the VtkVisImageItem class.
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

// ** INCLUDES **
#include "VtkVisPipelineItem.h"

class vtkAlgorithm;
class vtkImageChangeInformation;
class vtkPointSet;
class vtkProp3D;
class vtkRenderer;

class VtkAlgorithmProperties;
class VtkCompositeFilter;

/**
 * \brief An item in the VtkVisPipeline containing an image to be visualized.
 *
 * Any vtkImageAlgorithm object is represented by a VtkVisImageItem to be assigned a mapper,
 * an actor and its visualization properties.
 * \sa VtkVisPipelineItem
 */
class VtkVisImageItem : public VtkVisPipelineItem
{
public:
    /// @brief Constructor for a source/filter object.
    VtkVisImageItem(vtkAlgorithm* algorithm,
                    TreeItem* parentItem,
                    const QList<QVariant> data = QList<QVariant>());

    /// @brief Constructor for composite filter
    VtkVisImageItem(VtkCompositeFilter* compositeFilter, TreeItem* parentItem,
                    const QList<QVariant> data = QList<QVariant>());

    ~VtkVisImageItem() override;

    /// @brief Initializes vtkMapper and vtkActor necessary for visualization of
    /// the item and sets the item's properties.
    void Initialize(vtkRenderer* renderer) override;

    void setTranslation(double x, double y, double z) const override;

    vtkAlgorithm* transformFilter() const override;

protected:
    /// Selects the appropriate VTK-Writer object and writes the object to a file with the given name.
    int callVTKWriter(vtkAlgorithm* algorithm,
                      const std::string& filename) const override;
    void setVtkProperties(VtkAlgorithmProperties* vtkProps);

private:
    vtkImageChangeInformation* _transformFilter;
};
