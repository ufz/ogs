/**
 * \file SurfaceModel.cpp
 * 23/04/2010 KR Initial implementation
 *
 */

#include "SurfaceModel.h"

#include "VtkSurfacesSource.h"

SurfaceModel::SurfaceModel( QString name, const GEOLIB::SurfaceVec* surfaceVec, QObject* parent /*= 0*/ )
: Model(name, parent), _surfaceVec(surfaceVec)
{
	QList<QVariant> rootData;
	delete _rootItem;
	rootData << "Id" << "x" << "y" << "z";
	_rootItem = new TreeItem(rootData, NULL);
	setData(_surfaceVec, _rootItem);

	this->constructVTKObject();
}

SurfaceModel::~SurfaceModel()
{
	_vtkSource->Delete();
}

int SurfaceModel::columnCount( const QModelIndex& parent /*= QModelIndex()*/ ) const
{
	Q_UNUSED(parent)

	return 4;
}

void SurfaceModel::constructVTKObject()
{
	_vtkSource = VtkSurfacesSource::New();
	VtkSurfacesSource* source = static_cast<VtkSurfacesSource*>(_vtkSource);
	source->SetName(this->_name + " - Surfaces");
	source->setSurfaces(_surfaceVec->getVector());
}

void SurfaceModel::setData(const GEOLIB::SurfaceVec* surfaceVec, TreeItem* parent)
{
	Q_UNUSED(parent)
	const std::vector<GEOLIB::Surface*> *surfaces = surfaceVec->getVector();
	int nSurfaces = surfaces->size();
	for (int i=0; i<nSurfaces; i++)
	{
		QList<QVariant> surface;
		std::string sfc_name("");
		surface << "Surface " + QString::number(i);
		if (surfaceVec->getNameOfElementByID(i, sfc_name)) surface << QString::fromStdString(sfc_name);

		TreeItem* surfaceItem = new TreeItem(surface, _rootItem);
		_rootItem->appendChild(surfaceItem);

		const std::vector<GEOLIB::Point*> *nodesVec = (*surfaces)[i]->getPointVec();

		int nElems = static_cast<int>((*surfaces)[i]->getNTriangles());
		for (int j=0; j<nElems; j++)
		{
			QList<QVariant> elem;
			elem << j << static_cast<int>((*(*(*surfaces)[i])[j])[0]) << static_cast<int>((*(*(*surfaces)[i])[j])[1]) << static_cast<int>((*(*(*surfaces)[i])[j])[2]);
			TreeItem* child = new TreeItem(elem, surfaceItem);
			surfaceItem->appendChild(child);

			for (int k=0; k<3; k++)
			{
				QList<QVariant> node;
				node << static_cast<int>((*(*(*surfaces)[i])[j])[k]) << QString::number((*(*nodesVec)[(*(*(*surfaces)[i])[j])[k]])[0],'f') << QString::number((*(*nodesVec)[(*(*(*surfaces)[i])[j])[k]])[1],'f') << QString::number((*(*nodesVec)[(*(*(*surfaces)[i])[j])[k]])[2],'f');
				TreeItem* nchild = new TreeItem(node, child);
				child->appendChild(nchild);
			}
		}
	}

	reset();
}

void SurfaceModel::updateData()
{
	clearData();
	this->_vtkSource->Delete();

	TreeModel::updateData();
	this->constructVTKObject();
}
