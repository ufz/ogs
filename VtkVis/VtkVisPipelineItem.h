/**
 * \file VtkVisPipelineItem.h
 * 17/2/2010 LB Initial implementation
 *
 */

#ifndef VTKVISPIPELINEITEM_H
#define VTKVISPIPELINEITEM_H

// ** INCLUDES **
#include "Configure.h"
#include "TreeItem.h"

#include <QList>
#include <QMap>
//#include <QObject>
#include <QString>
#include <QVariant>

class vtkAlgorithm;
class vtkPointSet;
class QVtkDataSetMapper;
class vtkProp3D;
class vtkRenderer;
class VtkAlgorithmProperties;
class vtkOsgActor;
class VtkCompositeFilter;
class vtkTransformFilter;
class vtkDataSetAttributes;

/**
 * \brief An item in the VtkVisPipeline containing a graphic object to be visualized.
 *
 * Any VTK-object (source-items, filter-items, etc.) need to be put into a VtkPipelineItem
 * to be assigned a mapper, an actor and its visualization properties (colour, etc.).
 */
class VtkVisPipelineItem : /*public QObject,*/ public TreeItem
{
//	Q_OBJECT

public:
	/// @brief Constructor for a source/filter object.
	VtkVisPipelineItem(vtkAlgorithm* algorithm,
	                   TreeItem* parentItem,
	                   const QList<QVariant> data = QList<QVariant>());

	/// @brief Constructor for composite filter
	VtkVisPipelineItem(VtkCompositeFilter* compositeFilter, TreeItem* parentItem,
	                   const QList<QVariant> data = QList<QVariant>());

	~VtkVisPipelineItem();

	/// @brief Returns a VtkVisPipelineItem.
	VtkVisPipelineItem* child(int row) const;

	/// @brief Initializes vtkMapper and vtkActor necessary for visualization of
	/// the item and sets the item's properties.
	virtual void Initialize(vtkRenderer* renderer) = 0;

	QVariant data(int column) const;
	bool setData(int column, const QVariant &value);

	/// @brief Returns the algorithm object
	vtkAlgorithm* algorithm() const { return _algorithm; }

	/// @brief Returns the actor as vtkProp3D
	vtkProp3D* actor() const;

	// Dummy for implementation in derived classes
	virtual const QString GetActiveAttribute() const { return QString(""); }

	// Dummy for implementation in derived classes
	virtual void SetActiveAttribute(const QString& str) { (void)str; }

	/// @brief Returns the mapper
	QVtkDataSetMapper* mapper() const { return _mapper; }

	/// @brief Returns the composite filter
	VtkCompositeFilter* compositeFilter() const { return _compositeFilter; }

	/// @brief Returns if the VTK object is visible in the visualization.
	bool isVisible() const;

	/// @brief Sets the visibility of the VTK object in the visualization.
	void setVisible(bool visible);

	/// @brief Writes this algorithm's vtkDataSet (i.e. vtkPolyData or vtkUnstructuredGrid) to a vtk-file.
	int writeToFile(const std::string &filename) const;

	/// @brief Sets the geometry and data scaling.
	virtual void setScale(double x, double y, double z) const;

	/// @brief Translates the item in vis-space.
	virtual void setTranslation(double x, double y, double z) const;

	// Dummy for implementation in derived classes
	virtual vtkTransformFilter* transformFilter() const { return NULL; }
	/// @brief Sets the geometry and date scaling recursively on all children of
	/// this item.
	void setScaleOnChildren(double x, double y, double z) const;

protected:
	vtkProp3D* _actor;
	vtkAlgorithm* _algorithm;
	QVtkDataSetMapper* _mapper;
	vtkRenderer* _renderer;
	VtkCompositeFilter* _compositeFilter;

	virtual int callVTKWriter(vtkAlgorithm* algorithm, const std::string &filename) const;

	void SetScalarVisibility(bool on);

private:
};

#endif // VTKVISPIPELINEITEM_H

