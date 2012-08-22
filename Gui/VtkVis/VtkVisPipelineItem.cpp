/**
 * \file VtkVisPipelineItem.cpp
 * 17/2/2010 LB Initial implementation
 *
 * Implementation of VtkVisPipelineItem
 */

// ** INCLUDES **
#include "VtkAlgorithmProperties.h"
#include "VtkVisPipelineItem.h"
#include "VtkCompositeFilter.h"

#include "QVtkDataSetMapper.h"
#include <vtkActor.h>
#include <vtkAlgorithm.h>
#include <vtkDataSetMapper.h>
#include <vtkPointSet.h>
#include <vtkProperty.h>
#include <vtkRenderer.h>
#include <vtkSmartPointer.h>
#include <vtkTextureMapToPlane.h>
#include <vtkGenericDataObjectWriter.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkImageActor.h>

#include <QMessageBox>

#ifdef OGS_USE_OPENSG
#include "vtkOsgConverter.h"
#include <OpenSG/OSGSceneFileHandler.h>
#endif

VtkVisPipelineItem::VtkVisPipelineItem(
        vtkAlgorithm* algorithm, TreeItem* parentItem,
        const QList<QVariant> data /*= QList<QVariant>()*/)
	: TreeItem(data, parentItem),   _actor(NULL), _algorithm(algorithm),
	  _renderer(NULL), _compositeFilter(NULL), _vtkProps(NULL)
{
	VtkVisPipelineItem* visParentItem = dynamic_cast<VtkVisPipelineItem*>(parentItem);
	if (parentItem->parentItem())
		_algorithm->SetInputConnection(visParentItem->algorithm()->GetOutputPort());
}

VtkVisPipelineItem::VtkVisPipelineItem(
        VtkCompositeFilter* compositeFilter, TreeItem* parentItem,
        const QList<QVariant> data /*= QList<QVariant>()*/)
	: TreeItem(data, parentItem), _actor(NULL), _renderer(NULL), _compositeFilter(compositeFilter),
	  _vtkProps(NULL)
{
	_algorithm = _compositeFilter->GetOutputAlgorithm();
}

VtkVisPipelineItem::~VtkVisPipelineItem()
{
	_renderer->RemoveActor(_actor);
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
		return isVisible();
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
			if(!dynamic_cast<vtkImageActor*>(_actor))
			{
				vtkOsgConverter osgConverter(static_cast<vtkActor*>(_actor));
				if(osgConverter.WriteAnActor())
					OSG::SceneFileHandler::the().write(
					        osgConverter.GetOsgNode(), filename.c_str());
			}
			else
				QMessageBox::warning(NULL, "Conversion to OpenSG not possible",
					"It is not possible to convert an vtkImageData based object\nto OpenSG. If you want to convert raster data import it via \" File / Import / Raster Files as PolyData\"!");
#else       	
			QMessageBox::warning(
				NULL,
				"Functionality not implemented",
				"Sorry but this program was not compiled with OpenSG support.");
#endif
			return 0;
		}

		return callVTKWriter(this->algorithm(), filename);
	}
	return 0;
}

int VtkVisPipelineItem::callVTKWriter(vtkAlgorithm* algorithm, const std::string &filename) const
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

void VtkVisPipelineItem::setScale(double x, double y, double z) const
{
	(void)x;
	(void)y, (void)z;
}

void VtkVisPipelineItem::setTranslation(double x, double y, double z) const
{
	(void)x;
	(void)y, (void)z;
}

void VtkVisPipelineItem::setScaleOnChildren(double x, double y, double z) const
{
	for (int i = 0; i < this->childCount(); ++i)
	{
		VtkVisPipelineItem* child = this->child(i);
		child->setScale(x, y, z);
	}
}

void VtkVisPipelineItem::setBackfaceCulling(bool enable) const
{
	// Reimplemented in subclass
	(void)enable;
}

void VtkVisPipelineItem::setBackfaceCullingOnChildren(bool enable) const
{
	for (int i = 0; i < this->childCount(); ++i)
	{
		VtkVisPipelineItem* child = this->child(i);
		child->setBackfaceCulling((int)enable);
		child->setBackfaceCullingOnChildren((int)enable);
	}
}

QStringList VtkVisPipelineItem::getScalarArrayNames() const
{
	vtkDataSet* dataSet = vtkDataSet::SafeDownCast(this->algorithm()->GetOutputDataObject(0));
	QStringList dataSetAttributesList;
	if (dataSet)
	{
		vtkPointData* pointData = dataSet->GetPointData();
		//std::cout << "  #point data arrays: " << pointData->GetNumberOfArrays() << std::endl;
		for (int i = 0; i < pointData->GetNumberOfArrays(); i++)
			//std::cout << "    Name: " << pointData->GetArrayName(i) << std::endl;
			dataSetAttributesList.push_back(QString("P-") + pointData->GetArrayName(i));

		vtkCellData* cellData = dataSet->GetCellData();
		//std::cout << "  #cell data arrays: " << cellData->GetNumberOfArrays() << std::endl;
		for (int i = 0; i < cellData->GetNumberOfArrays(); i++)
			//std::cout << "    Name: " << cellData->GetArrayName(i) << std::endl;
			dataSetAttributesList.push_back(QString("C-") + cellData->GetArrayName(i));
	}
	return dataSetAttributesList;
}