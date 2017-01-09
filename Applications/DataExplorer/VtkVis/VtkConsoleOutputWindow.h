/**
* \brief  VtkWin32ConsoleOutputWindow is used to suppress message boxes on Windows.
*
* \copyright
* Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
*            Distributed under a Modified BSD License.
*              See accompanying file LICENSE.txt or
*              http://www.opengeosys.org/project/license
*
*/

#ifndef VTKCONSOLEOUTPUTWINDOW_H_
#define VTKCONSOLEOUTPUTWINDOW_H_

#include "vtkOutputWindow.h"

/// Is used to suppress message boxes on Windows and instead print VTK warnings
/// and errors in the console. Suppresses some specific messages as defined in
/// VtkConsoleOutputWindow::DisplayText()
class VtkConsoleOutputWindow : public vtkOutputWindow
{
public:
    vtkTypeMacro(VtkConsoleOutputWindow,vtkOutputWindow);
    void PrintSelf(ostream& os, vtkIndent indent) override;

    static VtkConsoleOutputWindow * New();
    virtual void DisplayText(const char*);

protected:
    VtkConsoleOutputWindow();
    virtual ~VtkConsoleOutputWindow();

private:
    VtkConsoleOutputWindow(const VtkConsoleOutputWindow &);  // Not implemented.
    void operator=(const VtkConsoleOutputWindow &);  // Not implemented.
};

#endif // VTKCONSOLEOUTPUTWINDOW_H_
