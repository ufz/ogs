/**
 * \file pymainwindow.cpp
 * 2/6/2010 LB Initial implementation
 * 
 * Implementation of pymainwindow
 */

// ** INCLUDES **
#include "boost/python.hpp"

#include "mainwindow.h"

using namespace boost::python;
namespace python = boost::python;

BOOST_PYTHON_MODULE(ogsgui)
{
	class_<MainWindow, boost::noncopyable>("ogsgui", init<>())
		.def("ShowWindow", &MainWindow::ShowWindow)
		.def("HideWindow", &MainWindow::HideWindow)
	;

	class_<StartQt4, boost::noncopyable>("StartQt4", init<>());
}
/**
TODO
vtkWidget entfernen, sonst Konflikt zwischen VRED Renderfenster und vtkWidget

Weiteres Problem: vred Prozess läuft nach Beenden von VRED weiter

import imp
ogsmodule = imp.load_dynamic('ogsgui', 'E:/bilke/geosys/branch/sources/Build/lib/Release/ogs-gui-vred.dll')
qt = ogsmodule.StartQt4()
ogs = ogsmodule.ogsgui()
ogs.ShowWindow()
 */