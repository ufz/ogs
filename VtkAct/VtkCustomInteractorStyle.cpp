/**
 * \file VtkInteractorStyle.cpp
 * 21/6/2010 LB Initial implementation
 * 
 * Implementation of VtkInteractorStyle
 */

// ** INCLUDES **
#include "VtkCustomInteractorStyle.h"

#include <vtkRenderWindowInteractor.h>
#include <vtkObjectFactory.h>
#include <vtkProp.h>

#include <string>

vtkStandardNewMacro(VtkCustomInteractorStyle);

VtkCustomInteractorStyle::VtkCustomInteractorStyle()
: _highlightActor(true)
{

}

void VtkCustomInteractorStyle::OnChar()
{
	switch (Interactor->GetKeyCode())
	{
	case '3':
		std::cout << "The 3 key was pressed." << std::endl;
		break;
	case 'a':
		break;
	default:
		vtkInteractorStyleTrackballCamera::OnChar();
	}
}

void VtkCustomInteractorStyle::highlightActor( vtkProp3D* actor )
{
	if (_highlightActor)
		HighlightProp((vtkProp*)actor);
}

void VtkCustomInteractorStyle::setHighlightActor(bool on)
{ 
	_highlightActor = on;
	if (!on)
		HighlightProp((vtkProp*)NULL);
}
