/**
 * \file
 * \author Lars Bilke
 * \date   2010-10-25
 * \brief  Implementation of the vtk to osg converter tool.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

// ** INCLUDES **
#include "VtkOsgConverter.h"

#include <iostream>

#include <vtkActor.h>
#include <vtkDataSetMapper.h>
#include <vtkGenericDataObjectReader.h>
#include <vtkGeometryFilter.h>
#include <vtkImageDataGeometryFilter.h>
#include <vtkImageMapper.h>
#include <vtkPolyDataMapper.h>
#include <vtkPolyDataNormals.h>
#include <vtkSmartPointer.h>
#include <vtkXMLImageDataReader.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkXMLRectilinearGridReader.h>
#include <vtkXMLStructuredGridReader.h>
#include <vtkXMLUnstructuredGridReader.h>

#include <OpenSG/OSGSceneFileHandler.h>
#include <OpenSG/OSGSwitch.h>

#include <boost/filesystem/operations.hpp>
#include <boost/regex.hpp>
#include <vector>
using namespace boost::filesystem;
using namespace std;

// Replace file extension
void replaceExt(string& s, const string& newExt)
{
	string::size_type i = s.rfind('.', s.length());
	if (i != string::npos)
		s.replace(i + 1, newExt.length(), newExt);
}

// Get file extension
string getFileExt(const string& s)
{
	size_t i = s.rfind('.', s.length());
	if (i != string::npos)
		return s.substr(i + 1, s.length() - i);
	return "";
}

// No arguments: batch convert all vt* files
// switch argument: batch convert all vt* files into one osb file with a switch
// file argument: convert only the specified file
int main (int argc, char const* argv[])
{
	vector<string> filenames;
	bool useSwitch = false;
	if (argc == 2)
	{
		if (string(argv[1]).find("switch") != string::npos)
			useSwitch = true;
		else
			filenames.push_back(string(argv[1]));
	}

	if (useSwitch || filenames.empty())
	{
		const boost::regex e(".+\\.vt[a-z]");
		directory_iterator end;
		for (directory_iterator it("./"); it != end; ++it)
		{
			string curFile = it->path().filename().string();
			if (regex_match(curFile, e))
				filenames.push_back(curFile);
		}
	}

	OSG::osgInit(0, NULL);

	vtkPolyDataMapper* mapper = vtkPolyDataMapper::New();
	OSG::NodePtr switchNode = OSG::Node::create();
	OSG::SwitchPtr switchCore = OSG::Switch::create();
	beginEditCP(switchCore);
	switchCore->setChoice(0);
	endEditCP(switchCore);
	beginEditCP(switchNode);
	switchNode->setCore(switchCore);
	endEditCP(switchNode);

	for (vector<string>::const_iterator it = filenames.begin(); it != filenames.end(); ++it)
	{
		string filename(*it);
		cout << "Opening file " << filename << " ... " << endl << flush;
		string fileExt = getFileExt(filename);

		vtkXMLDataReader* reader = NULL;
		vtkGenericDataObjectReader* oldStyleReader = NULL;
		if (fileExt.find("vti") != string::npos)
		{
			reader = vtkXMLImageDataReader::New();
			vtkSmartPointer<vtkImageDataGeometryFilter> geoFilter =
			        vtkSmartPointer<vtkImageDataGeometryFilter>::New();
			geoFilter->SetInputConnection(reader->GetOutputPort());
			mapper->SetInputConnection(geoFilter->GetOutputPort());
		}
		if (fileExt.find("vtr") != string::npos)
		{
			reader = vtkXMLRectilinearGridReader::New();
			vtkSmartPointer<vtkGeometryFilter> geoFilter =
			        vtkSmartPointer<vtkGeometryFilter>::New();
			geoFilter->SetInputConnection(reader->GetOutputPort());
			mapper->SetInputConnection(geoFilter->GetOutputPort());
		}
		else if (fileExt.find("vts") != string::npos)
		{
			reader = vtkXMLStructuredGridReader::New();
			vtkSmartPointer<vtkGeometryFilter> geoFilter =
			        vtkSmartPointer<vtkGeometryFilter>::New();
			geoFilter->SetInputConnection(reader->GetOutputPort());
			mapper->SetInputConnection(geoFilter->GetOutputPort());
		}
		else if (fileExt.find("vtp") != string::npos)
		{
			reader = vtkXMLPolyDataReader::New();
			mapper->SetInputConnection(reader->GetOutputPort());
		}
		else if (fileExt.find("vtu") != string::npos)
		{
			reader = vtkXMLUnstructuredGridReader::New();
			vtkSmartPointer<vtkGeometryFilter> geoFilter =
			        vtkSmartPointer<vtkGeometryFilter>::New();
			geoFilter->SetInputConnection(reader->GetOutputPort());
			mapper->SetInputConnection(geoFilter->GetOutputPort());
		}
		else if (fileExt.find("vtk") != string::npos)
		{
			oldStyleReader = vtkGenericDataObjectReader::New();
			oldStyleReader->SetFileName(filename.c_str());
			oldStyleReader->Update();
			if(oldStyleReader->IsFilePolyData())
				mapper->SetInputConnection(oldStyleReader->GetOutputPort());
			else
			{
				vtkSmartPointer<vtkGeometryFilter> geoFilter =
				        vtkSmartPointer<vtkGeometryFilter>::New();
				geoFilter->SetInputConnection(oldStyleReader->GetOutputPort());
				mapper->SetInputConnection(geoFilter->GetOutputPort());
			}
		}
		else
		{
			cout << "Not a valid vtk file ending (vti, vtr, vts, vtp, vtu, vtk)" <<
			endl;
			return 1;
		}

		if (fileExt.find("vtk") == string::npos)
		{
			reader->SetFileName(filename.c_str());
			reader->Update();
		}

		vtkActor* actor = vtkActor::New();
		actor->SetMapper(mapper);

		vtkOsgConverter converter(actor);
		converter.SetVerbose(true);
		//converter->SetMapper(mapper);
		converter.WriteAnActor();
		OSG::NodePtr node = converter.GetOsgNode();
		replaceExt(filename, "osb");
		if (useSwitch)
		{
			beginEditCP(switchNode);
			switchNode->addChild(node);
			endEditCP(switchNode);
		}
		else
			OSG::SceneFileHandler::the().write(node, filename.c_str());

		if (reader)
			reader->Delete();
		if (oldStyleReader)
			oldStyleReader->Delete();
	}
	if (useSwitch)
	{
		string filename(filenames[0]);
		replaceExt(filename, "osb");
		OSG::SceneFileHandler::the().write(switchNode, filename.c_str());
	}
	//mapper->Delete(); // TODO crashes

	OSG::osgExit();

	cout << "File conversion finished" << endl;

	return 0;
}
