/**
 * \file VtkVisPipelineItem.cpp
 * 17/2/2010 LB Initial implementation
 *
 * Implementation of VtkVisPipelineItem
 */

// ** INCLUDES **
#include "VtkVisPipelineItem.h"
#include "VtkAlgorithmProperties.h"

#include <vtkAlgorithm.h>
#include <vtkPointSet.h>
#include <vtkDataSetMapper.h>
#include "QVtkDataSetMapper.h"
#include <vtkActor.h>
#include <vtkRenderer.h>
#include <vtkProperty.h>
#include <vtkSmartPointer.h>
#include <vtkTextureMapToPlane.h>

#include <vtkGenericDataObjectWriter.h>

#include <QMessageBox>

#include "VtkCompositeFilter.h"

#include <vtkPointData.h>
#include <vtkCellData.h>

#ifdef OGS_USE_OPENSG
#include "vtkOsgConverter.h"
#include <OpenSG/OSGSceneFileHandler.h>
#endif


VtkVisPipelineItem::VtkVisPipelineItem(
	vtkAlgorithm* algorithm, TreeItem* parentItem,
	const QList<QVariant> data /*= QList<QVariant>()*/)
: TreeItem(data, parentItem),	_actor(NULL), _algorithm(algorithm), _mapper(NULL), _renderer(NULL),
	  _compositeFilter(NULL)
{
	VtkVisPipelineItem* visParentItem = dynamic_cast<VtkVisPipelineItem*>(parentItem);
	if (parentItem->parentItem())
	{
		_algorithm->SetInputConnection(visParentItem->algorithm()->GetOutputPort());
	}
}

VtkVisPipelineItem::VtkVisPipelineItem(
	VtkCompositeFilter* compositeFilter, TreeItem* parentItem,
	const QList<QVariant> data /*= QList<QVariant>()*/)
: TreeItem(data, parentItem), 	_actor(NULL), _mapper(NULL), _renderer(NULL),
	  _compositeFilter(compositeFilter)
{
	_algorithm = _compositeFilter->GetOutputAlgorithm();
}

VtkVisPipelineItem::~VtkVisPipelineItem()
{
	_renderer->RemoveActor(_actor);
	_mapper->Delete();
	_actor->Delete();
	delete _compositeFilter;
}

VtkVisPipelineItem* VtkVisPipelineItem::child( int row ) const
{
	TreeItem* treeItem = TreeItem::child(row);
	if (treeItem)
		return dynamic_cast<VtkVisPipelineItem*>(treeItem);
	else
		return NULL;
}

QVariant VtkVisPipelineItem::data( int column ) const
{
	if (column == 1)
	{
		return isVisible();
	}
	else
		return TreeItem::data(column);
}

bool VtkVisPipelineItem::setData( int column, const QVariant &value )
{
	if (column == 1)
	{
		setVisible(value.toBool());
		return true;
	}
	else
		return TreeItem::setData(column, value);

}
bool VtkVisPipelineItem::isVisible() const
{
	return (bool)_actor->GetVisibility();
}

void VtkVisPipelineItem::setVisible( bool visible )
{
	_actor->SetVisibility((int)visible);
	_actor->Modified();
	_renderer->Render();
}

int VtkVisPipelineItem::writeToFile(const std::string &filename) const
{
	if (!filename.empty())
	{
		if (filename.substr(filename.size() - 4).find("os") != std::string::npos)
		{
			#ifdef OGS_USE_OPENSG
			vtkOsgConverter osgConverter(static_cast<vtkActor*>(_actor));
			if(osgConverter.WriteAnActor())
				OSG::SceneFileHandler::the().write(osgConverter.GetOsgNode(), filename.c_str());
			#else
			QMessageBox::warning(NULL, "Functionality not implemented",
				"Sorry but this program was not compiled with OpenSG support.");
			#endif
			return 0;
		}

		return callVTKWriter(this->algorithm(), filename);
	}
	return 0;
}

int VtkVisPipelineItem::callVTKWriter(const vtkAlgorithm* algorithm, const std::string &filename) const
{
	// needs to be implemented in derived classes!
	(void)algorithm;
	(void)filename;
	return 0;
}

vtkProp3D* VtkVisPipelineItem::actor() const
{
	return _actor;
}

void VtkVisPipelineItem::SetScalarVisibility( bool on )
{
	_mapper->SetScalarVisibility(on);
}

void VtkVisPipelineItem::setScale(double x, double y, double z) const
{
	(void)x; (void)y, (void)z;
}

void VtkVisPipelineItem::setTranslation(double x, double y, double z) const
{
	(void)x; (void)y, (void)z;
}

void VtkVisPipelineItem::setScaleOnChildren(double x, double y, double z) const
{
	for (int i = 0; i < this->childCount(); ++i)
	{
		VtkVisPipelineItem* child = this->child(i);
		child->setScale(x, y, z);
	}
}
