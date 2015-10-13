/**
*
* \copyright
* Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
*            Distributed under a Modified BSD License.
*              See accompanying file LICENSE.txt or
*              http://www.opengeosys.org/project/license
*
*/
#include "VtkConsoleOutputWindow.h"

#include <string>

#include "vtkObjectFactory.h"
#ifdef WIN32
#include "vtkWindows.h"
#endif

vtkStandardNewMacro(VtkConsoleOutputWindow);

//----------------------------------------------------------------------------
VtkConsoleOutputWindow::VtkConsoleOutputWindow()
{

}

//----------------------------------------------------------------------------
VtkConsoleOutputWindow::~VtkConsoleOutputWindow()
{
}

//----------------------------------------------------------------------------
// Display text in the window, and translate the \n to \r\n.
//
void VtkConsoleOutputWindow::DisplayText(const char* someText)
{
	if(!someText)
		return;

	// Disable warnings
	std::string someTextString(someText);
	if((someTextString.find(
		"This is very expensive for vtkMappedDataArray subclasses, since the scalar array must be generated for each call.") != std::string::npos) ||
	   (someTextString.find("Invalid framebuffer operation") != std::string::npos))
		return;

	// Create a buffer big enough to hold the entire text
	char* buffer = new char[strlen(someText)+1];
	// Start at the beginning
	const char* NewLinePos = someText;
	while(NewLinePos)
	{
		int len = 0;
		// Find the next new line in text
		NewLinePos = strchr(someText, '\n');
		// if no new line is found then just add the text
		if(NewLinePos == 0)
		{
#ifdef WIN32
			OutputDebugString(someText);
#endif
			cerr << someText;
		}
			// if a new line is found copy it to the buffer
			// and add the buffer with a control new line
		else
		{
			len = NewLinePos - someText;
			strncpy(buffer, someText, len);
			buffer[len] = 0;
			someText = NewLinePos+1;
#ifdef WIN32
			OutputDebugString(buffer);
			OutputDebugString("\r\n");
#endif
			cerr << buffer;
			cerr << "\r\n";
		}
	}
	delete [] buffer;
}

//----------------------------------------------------------------------------
void VtkConsoleOutputWindow::PrintSelf(ostream& os, vtkIndent indent)
{
	this->Superclass::PrintSelf(os,indent);
}
