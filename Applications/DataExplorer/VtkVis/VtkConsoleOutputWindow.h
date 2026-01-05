// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

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
    void DisplayText(const char* /*unused*/) override;

    VtkConsoleOutputWindow(const VtkConsoleOutputWindow&) = delete;
    void operator=(const VtkConsoleOutputWindow&) = delete;

protected:
    VtkConsoleOutputWindow();
    ~VtkConsoleOutputWindow() override;
};
