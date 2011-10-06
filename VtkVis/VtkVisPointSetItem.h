/**
 * \file VtkVisPointSetItem.h
 * 2011/09/29 KR Initial implementation
 *
 */


#ifndef VTKVISPOINTSETITEM_H
#define VTKVISPOINTSETITEM_H

// ** INCLUDES **
#include "VtkVisPipelineItem.h"

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
class VtkVisPointSetItem : public VtkVisPipelineItem
{
//	Q_OBJECT

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
	const QString GetActiveAttribute() const {return _activeAttribute; }

	/// @brief Initializes vtkMapper and vtkActor necessary for visualization of
	/// the item and sets the item's properties.
	void Initialize(vtkRenderer* renderer);

	vtkTransformFilter* transformFilter() const { return _transformFilter; }

	/// @brief Sets the selected attribute array for the visualisation of the data set.
	void SetActiveAttribute(const QString& name);

	void SetScalarRange(double min, double max);

	/// @brief Sets the geometry and data scaling.
	void setScale(double x, double y, double z) const;

	/// @brief Translates the item in vis-space.
	void setTranslation(double x, double y, double z) const;

protected:
	vtkTransformFilter* _transformFilter;
	QString _activeAttribute;

	virtual int callVTKWriter(vtkAlgorithm* algorithm, const std::string &filename) const;

	/// Sets a color lookup table for the current scalar array.
	void setLookupTableForActiveScalar();

	/// @brief Sets pre-set properties on vtkActor and on vtkMapper
	void setVtkProperties(VtkAlgorithmProperties* vtkProps);


private:
	/// @see SetActiveAttribute()
	bool setActiveAttributeOnData(vtkDataSetAttributes* data, std::string& name);

};

#endif // VTKVISPOINTSETITEM_H

