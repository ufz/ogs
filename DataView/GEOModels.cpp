/**
 * \file GEOModels.cpp
 * 9/2/2010 LB Initial implementation
 *
 * Implementation of GEOModels
 */

// ** INCLUDES **
#include "GEOModels.h"

#include "Model.h"
#include "PntsModel.h"
#include "StationTreeModel.h"
#include "LinesModel.h"
#include "SurfaceModel.h"

GEOModels::GEOModels( QObject* parent /*= 0*/ ) :
	QObject (parent)
{
	_stationModel = new StationTreeModel();
}

GEOModels::~GEOModels()
{
	delete _stationModel;
}

void GEOModels::addPointVec( std::vector<GEOLIB::Point*> *points, std::string &name, std::map<std::string, size_t>* name_pnt_id_map )
{
	GEOObjects::addPointVec(points, name, name_pnt_id_map);

	PntsModel* model = new PntsModel(QString::fromStdString(name), this->getPointVecObj(name), this);
	_pntModels.push_back(model);
	emit pointModelAdded(model);
}

bool GEOModels::appendPointVec(const std::vector<GEOLIB::Point*> &points, std::string &name)
{
	bool ret (GEOLIB::GEOObjects::appendPointVec (points, name));

	// search model
	QString qname (name.c_str());
	bool nfound (true);
	std::vector<PntsModel*>::iterator it(_pntModels.begin());
	while (nfound && it != _pntModels.end()) {
		if (((*it)->name()).contains (qname)) nfound = false;
		else it++;
	}
	if (nfound) std::cerr << "Model not found" << std::endl;
	else (*it)->updateData ();

	return ret;
}

bool GEOModels::removePointVec( const std::string &name )
{
	if (! isPntVecUsed(name)) {
		for (std::vector<PntsModel*>::iterator it = _pntModels.begin(); it
				!= _pntModels.end(); ++it) {
			if ((*it)->name().toStdString() == name) {
				emit pointModelRemoved(*it);
				delete *it;
				_pntModels.erase(it);
				break;
			}
		}
		return GEOObjects::removePointVec(name);
	}
	std::cout << "GEOModels::removePointVec() - There are still Polylines or Surfaces depending on these points." << std::endl;
	return false;
}

void GEOModels::addStationVec( std::vector<GEOLIB::Point*> *stations, std::string &name, const GEOLIB::Color* const color )
{
	GEOObjects::addStationVec(stations, name, color);

	_stationModel->addStationList(QString::fromStdString(name), stations);
	emit stationVectorAdded(_stationModel, name);
}

void GEOModels::filterStationVec(const std::string &name, const std::vector<PropertyBounds> &bounds)
{
	emit stationVectorRemoved(_stationModel, name);
	const std::vector<GEOLIB::Point*> *stations (GEOObjects::getStationVec(name));
	_stationModel->filterStations(name, stations, bounds);
	emit stationVectorAdded(_stationModel, name);
}

bool GEOModels::removeStationVec( const std::string &name )
{
	emit stationVectorRemoved(_stationModel, name);
	_stationModel->removeStationList(name);
	return GEOObjects::removeStationVec(name);
}

void GEOModels::addPolylineVec( std::vector<GEOLIB::Polyline*> *lines, const std::string &name, std::map<std::string,size_t>* ply_names )
{
	GEOObjects::addPolylineVec(lines, name, ply_names);
	if (lines->empty()) return;

	PolylinesModel* model = new PolylinesModel(QString::fromStdString(name), this->getPolylineVecObj(name), this);
	_lineModels.push_back(model);

	connect(model, SIGNAL(requestAppendPolylines(const std::vector<GEOLIB::Polyline*>&, std::string&)),
		this, SLOT(appendPolylineVec(const std::vector<GEOLIB::Polyline*>&, std::string&)));

	emit polylineModelAdded(model);
}

bool GEOModels::appendPolylineVec(const std::vector<GEOLIB::Polyline*> &polylines, std::string &name)
{
	bool ret (GEOLIB::GEOObjects::appendPolylineVec (polylines, name));

	// search model
	QString qname (name.c_str());
	bool nfound (true);
	std::vector<PolylinesModel*>::iterator it(_lineModels.begin());
	while (nfound && it != _lineModels.end()) {
		if (((*it)->name()).contains (qname)) nfound = false;
		else it++;
	}
	if (nfound) std::cerr << "Model not found" << std::endl;
	else (*it)->updateData ();

	return ret;
}

bool GEOModels::removePolylineVec( const std::string &name )
{
	for (std::vector<PolylinesModel*>::iterator it = _lineModels.begin();
			it != _lineModels.end(); ++it) {
		if ((*it)->name().toStdString() == name) {
			emit polylineModelRemoved(*it);
			delete *it;
			_lineModels.erase(it);
			return GEOObjects::removePolylineVec (name);
		}
	}
	return false;
}

void GEOModels::addSurfaceVec( std::vector<GEOLIB::Surface*> *surfaces, const std::string &name, std::map<std::string,size_t>* sfc_names )
{
	GEOObjects::addSurfaceVec(surfaces, name, sfc_names);
	if (surfaces->empty()) return;

	SurfaceModel* model = new SurfaceModel(QString::fromStdString(name), this->getSurfaceVecObj(name), this);
	_surfaceModels.push_back(model);
	emit surfaceModelAdded(model);
}

bool GEOModels::removeSurfaceVec( const std::string &name )
{
	for (std::vector<SurfaceModel*>::iterator it = _surfaceModels.begin();
			it != _surfaceModels.end(); ++it) {
		if ((*it)->name().toStdString() == name) {
			emit surfaceModelRemoved(*it);
			delete *it;
			_surfaceModels.erase(it);
			return GEOObjects::removeSurfaceVec (name);
		}
	}
	return false;
}
