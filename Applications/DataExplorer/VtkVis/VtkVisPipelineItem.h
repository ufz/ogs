/**
 * \file
 * \author Lars Bilke
 * \date   2010-02-17
 * \brief  Definition of the VtkVisPipelineItem class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <QList>
#include <QMap>
#include <QString>
#include <QVariant>

#include "BaseLib/BuildInfo.h"

#include "TreeItem.h"

class QStringList;
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
 * \brief An item in the VtkVisPipeline containing a graphic object to be visualized.
 *
 * Any VTK-object (source-items, filter-items, etc.) need to be put into a VtkPipelineItem
 * to be assigned a mapper, an actor and its visualization properties (colour, etc.).
 */
class VtkVisPipelineItem : /*public QObject,*/ public TreeItem
{
//    Q_OBJECT

public:
    /// @brief Constructor for a source/filter object.
    VtkVisPipelineItem(vtkAlgorithm* algorithm,
                       TreeItem* parentItem,
                       const QList<QVariant> data = QList<QVariant>());

    /// @brief Constructor for composite filter
    VtkVisPipelineItem(VtkCompositeFilter* compositeFilter, TreeItem* parentItem,
                       const QList<QVariant> data = QList<QVariant>());

    ~VtkVisPipelineItem() override;

    /// @brief Returns a VtkVisPipelineItem.
    VtkVisPipelineItem* child(int row) const;

    /// @brief Initializes vtkMapper and vtkActor necessary for visualization of
    /// the item and sets the item's properties.
    virtual void Initialize(vtkRenderer* renderer) = 0;

    QVariant data(int column) const override;
    bool setData(int column, const QVariant& value) override;

    /// @brief Returns the algorithm object
    vtkAlgorithm* algorithm() const { return _algorithm; }

    /// @brief Returns the actor as vtkProp3D
    vtkProp3D* actor() const;

    // Dummy for implementation in derived classes
    virtual const QString GetActiveAttribute() const { return QString(""); }

    // Dummy for implementation in derived classes
    virtual void SetActiveAttribute(const QString& str) { (void)str; }

    /// @brief Returns the composite filter
    VtkCompositeFilter* compositeFilter() const { return _compositeFilter; }

    /// @brief Returns if the VTK object is visible in the visualization.
    bool isVisible() const;

    /// @brief Sets the visibility of the VTK object in the visualization.
    void setVisible(bool visible);

    /// @brief Writes this algorithm's vtkDataSet (i.e. vtkPolyData or vtkUnstructuredGrid)
    /// to a vtk-file.
    int writeToFile(const std::string &filename) const;

    /**
      * @brief Scales the data in visualisation-space.
     * This function is empty and needs to be implemented by derived classes.
     */
    virtual void setScale(double x, double y, double z) const;

    /// @brief Sets the geometry and date scaling recursively on all children of
    /// this item.
    void setScaleOnChildren(double x, double y, double z) const;

    /**
      * @brief Translates the item in visualisation-space.
     * This function is empty and needs to be implemented by derived classes.
     */
    virtual void setTranslation(double x, double y, double z) const;

    /**
      * Returns the transform filter for the object.
     * This function needs to be implemented by derived classes.
     */
    virtual vtkAlgorithm* transformFilter() const = 0;

    /// @brief Enables / disables backface culling.
    virtual void setBackfaceCulling(bool enable) const;

    /// @brief Enables / disables backface culling on all children.
    void setBackfaceCullingOnChildren(bool enable) const;

    /// @brief Returns a list of array names prefixed with P- or C-
    /// for point and cell data.
    QStringList getScalarArrayNames() const;

    ///    @brief Returns the VtkAlgorithmProperties.
    VtkAlgorithmProperties* getVtkProperties() const { return _vtkProps; };

protected:
    vtkProp3D* _actor;
    vtkAlgorithm* _algorithm;
    vtkRenderer* _renderer;
    VtkCompositeFilter* _compositeFilter;

    /// @brief The active VtkAlgorithmProperties.
    /// From algorithm, compositeFilter, or copied from parent
    VtkAlgorithmProperties* _vtkProps;

    /**
     * Selects the appropriate VTK-Writer object and writes the object to a file with the given name.
     * This function is empty and needs to be implemented by derived classes.
     */
    virtual int callVTKWriter(vtkAlgorithm* algorithm, const std::string &filename) const;

private:
};
