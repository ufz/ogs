/**
 * \file QSpaceNavigatorClient.cpp
 * 31/08/2010 LB Initial implementation
 * 
 * Implementation of QSpaceNavigatorClient class
 */

// ** INCLUDES **
#include "QSpaceNavigatorClient.h"

QSpaceNavigatorClient* QSpaceNavigatorClient::_singleton = 0;

QSpaceNavigatorClient::QSpaceNavigatorClient(QObject* parent /*= NULL*/)
: QObject(parent), SpaceNavigatorClient()
{
	//SpaceNavigatorClient::Instance();
	_time.start();
	_timer = new QTimer(this);
	connect( _timer, SIGNAL(timeout()), this, SLOT(update()) );
}

QSpaceNavigatorClient::~QSpaceNavigatorClient()
{
	
}

QSpaceNavigatorClient* QSpaceNavigatorClient::Instance(QObject* parent /*= NULL*/)
{
	if(_singleton == 0)
		_singleton = new QSpaceNavigatorClient(parent);
	return _singleton;
}

void QSpaceNavigatorClient::init(const char *deviceName,
	int updateInterval /*= 100*/, SpaceNavigatorAxes axis /*= Y*/)
{
	SpaceNavigatorClient::init(deviceName, axis);
	if (updateInterval > 0)
		_timer->start(updateInterval);
}

void QSpaceNavigatorClient::update()
{
	SpaceNavigatorClient::update();
	if (_unconsumedData)
	{
		_unconsumedData = false;
		double x, y, z, rx, ry, rz;
		getTranslation(x, y, z);
		getRotation(rx, ry, rz);
		//emit updated(x, y, z, rx, ry, rz);
		emit translated(x, y, z);
	}
}

int QSpaceNavigatorClient::getFrameTime()
{
	return _time.restart();
}