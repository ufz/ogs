/**
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.net/LICENSE.txt
 *
 * \file VrpnClient.h
 *
 * Created on 2010-08-30 by Lars Bilke
 */

#ifndef VRPNCLIENT_H
#define VRPNCLIENT_H

#include <QObject>
#include <QVector>
#include <vrpn_Analog.h>

class QString;
class vrpn_Analog_Remote;

class VrpnClient : public QObject
{
	Q_OBJECT

public:
	VrpnClient(QString deviceName, int updateInterval = 100, QObject* parent = NULL);
	virtual ~VrpnClient();

	double getAnalog(int channel);

public slots:
	void update();

protected:
	QString _deviceName;
	vrpn_Analog_Remote* _vrpnAnalog;
	QVector<double> _analogData;

signals:
	void positionChanged(double x, double y, double z);
	void rotationChanged(double x, double y, double z);
};

#endif // VRPNCLIENT_H
