/**
 * \file VtkVisImageItem.h
 * 2011/09/29 KR Initial implementation
 *
 */


#ifndef VTKVISIMAGEITEM_H
#define VTKVISIMAGEITEM_H

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

/**
 * \brief An item in the VtkVisPipeline containing a graphic object to be visualized.
 *
 * Any VTK-object (source-items, filter-items, etc.) need to be put into a VtkPipelineItem
 * to be assigned a mapper, an actor and its visualization properties (colour, etc.).
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

	~VtkVisImageItem();

	/// @brief Initializes vtkMapper and vtkActor necessary for visualization of
	/// the item and sets the item's properties.
	void Initialize(vtkRenderer* renderer);

protected:
	virtual int callVTKWriter(vtkAlgorithm* algorithm, const std::string &filename) const;
	void setVtkProperties(VtkAlgorithmProperties* vtkProps);

private:

};

#endif // VTKVISIMAGEITEM_H

