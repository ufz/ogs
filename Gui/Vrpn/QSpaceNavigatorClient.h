/**
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.net/LICENSE.txt
 *
 * \file QSpaceNavigatorClient.h
 *
 * Created on 2010-08-31 by Lars Bilke
 */

#ifndef QSPACENAVIGATORCLIENT_H
#define QSPACENAVIGATORCLIENT_H

#include "SpaceNavigatorClient.h"
#include <QObject>
#include <QTime>
#include <QTimer>

/// @brief Is a qt specific implementation of SpaceNavigatorClient.
class QSpaceNavigatorClient : public QObject, public SpaceNavigatorClient
{
	Q_OBJECT

public:
	/// @brief Returns the singleton of this class
	static QSpaceNavigatorClient* Instance(QObject* parent = NULL);

	/// @brief Initializes the SpaceNavigator.
	/// Connects with the server and registers the callback handlers.
	/// It is possible to specify the z-axis as the up-axis.
	/// Default the y-axis is the up-axis.
	/// @param deviceName Example: "SpaceNav@viswork01.intern.ufz.de"
	/// @param updateInterval The update interval in ms. If 0 then you have to
	/// manually call update().
	/// @param axis The up axis.
	void init(const char* deviceName, int updateInterval = 100, SpaceNavigatorAxes axis = Y);

public slots:
	/// @brief Calls the SpaceNavigatorClient::update() method and emits updated()
	/// and translated() signals if there is new data.
	/// Is called automatically if an updateInterval was given in init().
	void update();

protected:
	/// @brief The constructor is protected because of the singleton
	/// design pattern.
	QSpaceNavigatorClient(QObject* parent = NULL);

	/// @brief Destructor.
	virtual ~QSpaceNavigatorClient();

	/// @brief Returns the elapsed time since the last function call in ms.
	int getFrameTime();

private:
	/// Used in getFrameTime() to compute the elapsed time since last call.
	QTime _time;

	/// Calls update(), see init().
	QTimer* _timer;

	/// @brief This one points to the class itself.
	/// You can use only one QSpaceNavigatorClient because it´s static.
	/// This is needed for the callback methods which only
	/// can access static members.
	static QSpaceNavigatorClient* _singleton;

signals:
	void updated(double x, double y, double z, double rx, double ry, double rz);

	/// Is emitted when there is new translation data.
	void translated(double x, double y, double z);
};

#endif // QSPACENAVIGATORCLIENT_H
