/**
*
* \copyright
* Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
*            Distributed under a Modified BSD License.
*              See accompanying file LICENSE.txt or
*              http://www.opengeosys.org/project/license
*
*/
#include "VtkConsoleOutputWindow.h"

#include <string>
#include <ostream>

#include "vtkObjectFactory.h"
#ifdef WIN32
#include "vtkWindows.h"
#endif

vtkStandardNewMacro(VtkConsoleOutputWindow);

//----------------------------------------------------------------------------
VtkConsoleOutputWindow::VtkConsoleOutputWindow() = default;

//----------------------------------------------------------------------------
VtkConsoleOutputWindow::~VtkConsoleOutputWindow() = default;

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

#ifdef WIN32
    OutputDebugString(someTextString.c_str());
#endif
    std::cerr << someText;
}

//----------------------------------------------------------------------------
void VtkConsoleOutputWindow::PrintSelf(ostream& os, vtkIndent indent)
{
    this->Superclass::PrintSelf(os,indent);
}
