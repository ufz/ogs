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

#ifndef QVRPNARTTRACKINGCLIENT_H
#define QVRPNARTTRACKINGCLIENT_H

#include "VrpnArtTrackingClient.h"
#include <QObject>

#include <QTimer>

class QString;

class QVrpnArtTrackingClient : public QObject, public VrpnArtTrackingClient
{
	Q_OBJECT

public:
	/// @brief Returns the singleton of this class
	static QVrpnArtTrackingClient* Instance(QObject* parent = NULL);

	/// @brief Initializes and starts the tracking.
	/// @param deviceName The name of the vrpn device @ vrpn server host name
	/// e.g. DTrack@visserv3.intern.ufz.de
	/// @param updateInterval The update interval in ms. If 0 then you have to
	/// manually call MainLoop().
	void StartTracking(const char* deviceName, int updateInterval = 100);

	/// @brief Returns the device name.
	QString deviceName() const { return _deviceName; }

	/// @ Returns the update interval.
	int updateInterval() const { return _updateInterval; }

public slots:
	/// @brief Calls the vrpn mainloop functions. Must be called once per frame.
	/// Is called automatically if an updateInterval was given in init().
	void MainLoop();

protected:
	QVrpnArtTrackingClient(QObject* parent = NULL);
	virtual ~QVrpnArtTrackingClient();

private:
	/// Calls MainLoop(), see init().
	QTimer* _timer;

	/// The device name, e.g. DTrack@141.65.34.36
	QString _deviceName;

	/// The update interval
	int _updateInterval;

	/// @brief This one points to the class itself.
	/// You can use only one QVrpnArtTrackingClient because it´s static.
	/// This is needed for the callback methods which only
	/// can access static members.
	static QVrpnArtTrackingClient* _singleton;

signals:
	void positionUpdated(double x, double y, double z);
};

#endif // QVRPNARTTRACKINGCLIENT_H
