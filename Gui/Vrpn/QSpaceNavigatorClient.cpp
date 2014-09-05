/**
 * \file
 * \author Lars Bilke
 * \date 2010-08-31
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.net/LICENSE.txt
 *
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

void QSpaceNavigatorClient::init(const char* deviceName,
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