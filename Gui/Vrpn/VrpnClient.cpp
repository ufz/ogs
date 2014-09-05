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
#include "VrpnClient.h"

#include <QTimer>
#include <vrpn_Analog.h>

void VRPN_CALLBACK handle_analog( void* userData, const vrpn_ANALOGCB a )
{
	int i;
	const char* name = (const char*)userData;

	printf("Analog %s:\n         %5.2f", name, a.channel[0]);
	for (i = 1; i < a.num_channel; i++)
		printf(", %5.2f", a.channel[i]);
	printf(" (%d chans)\n", a.num_channel);

	// if (a.num_channel >= 3)
	//              emit positionChanged(a.channel[1], a.channel[2], a.channel[3]);
	//      if (a.num_channel >= 6)
	//              emit rotationChanged(a.channel[4], a.channel[5], a.channel[6]);
}

VrpnClient::VrpnClient( QString deviceName, int updateInterval /*= 100*/, QObject* parent /*= NULL*/ )
	: QObject(parent)
{
	// Create vrpn analog device
	_deviceName = deviceName;
	_vrpnAnalog = new vrpn_Analog_Remote( deviceName.toStdString().c_str() );
	_vrpnAnalog->register_change_handler( 0, handle_analog );

	int numChannels = _vrpnAnalog->getNumChannels();
	_analogData.fill(0.0, numChannels);

	// Call the vrpn mainloop every updateInterval ms
	QTimer* vrpnTimer = new QTimer(this);
	connect(vrpnTimer, SIGNAL(timeout()), this, SLOT(update()));
	vrpnTimer->start(updateInterval);
}
VrpnClient::~VrpnClient()
{
	_vrpnAnalog->unregister_change_handler( 0, handle_analog );
	delete _vrpnAnalog;
}

double VrpnClient::getAnalog(int channel)
{
	return _analogData[channel];
}

void VrpnClient::update()
{
	_vrpnAnalog->mainloop();
}