/**
 * \file VtkVisPipelineItem.h
 * 17/2/2010 LB Initial implementation
 *
 */


#ifndef VTKVISPIPELINEITEM_H
#define VTKVISPIPELINEITEM_H

// ** INCLUDES **
#include "TreeItem.h"

#include <QList>
#include <QMap>
//#include <QObject>
#include <QString>
#include <QVariant>

#ifdef OGS_USE_OPENSG
	#include <OpenSG/OSGNode.h>
	#include <OpenSG/OSGRefPtr.h>
#endif // OGS_USE_OPENSG

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
	void Initialize(vtkRenderer* renderer);

	QVariant data(int column) const;
	bool setData(int column, const QVariant &value);

	/// @brief Returns the algorithm object
	vtkAlgorithm* algorithm() const { return _algorithm; }

	/// @brief Returns the actor as vtkProp3D
	vtkProp3D* actor() const;

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

	vtkTransformFilter* transformFilter() const { return _transformFilter; }

	/// @brief Sets the selected attribute array for the visualisation of the data set.
	void SetActiveAttribute(const QString& name);

	void SetScalarRange(double min, double max);

	/// @brief Gets the last selected attribute.
	const QString& GetActiveAttribute() const {return _activeAttribute; }

	/// @brief Sets the geometry and data scaling.
	void setScale(double x, double y, double z) const;

	/// @brief Sets the geometry and date scaling recursively on all children of
	/// this item.
	void setScaleOnChildren(double x, double y, double z) const;
	
#ifdef OGS_USE_OPENSG
	// HACK static rootNode is set by VtkVisPipeline constructor
	/// Do not use this variable except in VtkVisPipeline constructor!
	static OSG::NodePtr rootNode;

protected:
	vtkOsgActor* _actor;
	OSG::RefPtr<OSG::NodePtr> _parentNode;
#else // OGS_USE_OPENSG
protected:
	vtkProp3D* _actor;
#endif // OGS_USE_OPENSG
	vtkAlgorithm* _algorithm;
	QVtkDataSetMapper* _mapper;
	vtkRenderer* _renderer;
	VtkCompositeFilter* _compositeFilter;
	vtkTransformFilter* _transformFilter;
	QString _activeAttribute;

	/// Sets a color lookup table for the current scalar array.
	void setLookupTableForActiveScalar();

	/// @brief Sets pre-set properties on vtkActor and on vtkMapper
	void setVtkProperties(VtkAlgorithmProperties* vtkProps);

	void SetScalarVisibility(bool on);
	
private:
	/// @see SetActiveAttribute()
	bool setActiveAttributeOnData(vtkDataSetAttributes* data, std::string& name);

};

#endif // VTKVISPIPELINEITEM_H

