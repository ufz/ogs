// SpaceNavigatorClient.cpp
//
//
// Author: Lars Bilke

// ** INCLUDES **
#include "SpaceNavigatorClient.h"
#include <iostream>


SpaceNavigatorClient* SpaceNavigatorClient::_spacenavigator = 0;

SpaceNavigatorClient::SpaceNavigatorClient() 
{
	_button = NULL;
	_analog = NULL;
	_unconsumedData = false;
	_frameRotationFactor = 0.0;
	_frameTranslationFactor = 0.0;
	_upAxis = Y;

	// set all fields to start values
	_spacenavigator = this;
	for(int i = 0; i < 8; i++)
		buttons[0] = false;
	_defaultButtonBehaviour = true;
	x = y = z = rx = ry = rz = 0.0;
	_translationFactor = _rotationFactor = 1.0;

	_dominating = false;
	_invertAxes[0] = _invertAxes[2] = _invertAxes[3] = _invertAxes[5] = 1.0;
	_invertAxes[1] = _invertAxes[4] = -1.0;

	_mode = TRANSROT;
}

SpaceNavigatorClient::~SpaceNavigatorClient()
{
	// clean up
	if(_button)
	{
		_button->unregister_change_handler(NULL, (vrpn_BUTTONCHANGEHANDLER)SpaceNavigatorClient::_handleButtons);
		delete _button;
		_button = NULL;
	}
	if(_analog)
	{
		_analog->unregister_change_handler(NULL, (vrpn_ANALOGCHANGEHANDLER)SpaceNavigatorClient::_handleAnalogs);
		delete _analog;
		_analog = NULL;
	}
}

void SpaceNavigatorClient::init( const char *deviceName, SpaceNavigatorAxes axis /*= Y*/ )
{
	// Create new button_remote object that connects to the corresponding server object
	// at the specified computer and register the callback function
	if(_button)
	{
		_button->unregister_change_handler(NULL, (vrpn_BUTTONCHANGEHANDLER)SpaceNavigatorClient::_handleButtons);
		delete _button;
	}
	_button = new vrpn_Button_Remote(deviceName);
	_button->register_change_handler(NULL, (vrpn_BUTTONCHANGEHANDLER)SpaceNavigatorClient::_handleButtons);

	// Create new analog_remote object that connects to the corresponding server object
	// at the specified computer and register the callback function
	if(_analog)
	{
		_analog->unregister_change_handler(NULL, (vrpn_ANALOGCHANGEHANDLER)SpaceNavigatorClient::_handleAnalogs);
		delete _analog;
	}
	_analog = new vrpn_Analog_Remote(deviceName);
	_analog->register_change_handler(NULL, (vrpn_ANALOGCHANGEHANDLER)SpaceNavigatorClient::_handleAnalogs);

	// set up axis
	if(axis == Z)
		_upAxis = axis;
	else
		_upAxis = Y;
}

void SpaceNavigatorClient::getTranslation(double& retx, double& rety, double& retz) const
{
	retx = x * _frameTranslationFactor;
	rety = y * _frameTranslationFactor;
	retz = z * _frameTranslationFactor;
}

void SpaceNavigatorClient::getRotation(double& retx, double& rety, double& retz) const
{
	retx = rx * _frameTranslationFactor;
	rety = ry * _frameTranslationFactor;
	retz = rz * _frameTranslationFactor;
}

void SpaceNavigatorClient::update()
{
	double frameTime = getFrameTime(); // in ms
	_frameTranslationFactor =  _translationFactor * (frameTime / 1000.0);
	_frameRotationFactor = _rotationFactor * (frameTime / 800.0); 
	
	// process messages from server
	this->_button->mainloop();
	this->_analog->mainloop();
}

void SpaceNavigatorClient::setTranslationFactor(double factor)
{
	_translationFactor = factor;
}

void SpaceNavigatorClient::setRotationFactor(double factor)
{
	_rotationFactor = factor;
}

bool SpaceNavigatorClient::getZUpAxis()
{
	if (_upAxis == Z)
		return true;
	else
		return false;
}

void SpaceNavigatorClient::setZUpAxis( bool zUp )
{
	if (zUp)
		_upAxis = Z;
	else
		_upAxis = Y;
}
void SpaceNavigatorClient::setDomination(bool dominating)
{
	this->_dominating = dominating;

	#ifdef SPACENAVIGATOR_DEBUG_OUTPUT
		std::cout << "SpaceNavigator: Axis domination mode: ";
		if(dominating)
			std::cout << "ON" << std::endl;
		else
			std::cout << "OFF" << std::endl;
	#endif // SPACENAVIGATOR_DEBUG_OUTPUT
}

void SpaceNavigatorClient::switchDomination()
{
	this->_dominating = !this->_dominating;
	
	#ifdef SPACENAVIGATOR_DEBUG_OUTPUT
		std::cout << "SpaceNavigator: Axis domination mode: ";
		if(_dominating)
			std::cout << "ON" << std::endl;
		else
			std::cout << "OFF" << std::endl;
	#endif // SPACENAVIGATOR_DEBUG_OUTPUT
}

void SpaceNavigatorClient::setMode(SpaceNavigatorClient::SpaceNavigatorMode mode)
{
	this->_mode = mode;

	#ifdef SPACENAVIGATOR_DEBUG_OUTPUT
		std::cout << "SpaceNavigator: Transformation mode: ";
		if(mode == TRANSROT)
			std::cout << "TRANSLATION and ROTATION" << std::endl;
		else if(mode == TRANS)
			std::cout << "TRANSLATION only" << std::endl;
		else if(mode == ROT)
			std::cout << "ROTATION only" << std::endl;
	#endif // SPACENAVIGATOR_DEBUG_OUTPUT
}

void SpaceNavigatorClient::switchMode()
{
	this->_mode = (SpaceNavigatorClient::SpaceNavigatorMode)((this->_mode + 1) % 3);

	#ifdef SPACENAVIGATOR_DEBUG_OUTPUT
		std::cout << "SpaceNavigator: Transformation mode: ";
		if(_mode == TRANSROT)
			std::cout << "TRANSLATION and ROTATION" << std::endl;
		else if(_mode == TRANS)
			std::cout << "TRANSLATION only" << std::endl;
		else if(_mode == ROT)
			std::cout << "ROTATION only" << std::endl;
	#endif // SPACENAVIGATOR_DEBUG_OUTPUT
}

void SpaceNavigatorClient::setDefaultButtonBehaviour(bool enable)
{
	_defaultButtonBehaviour = enable;
}


void SpaceNavigatorClient::invertAxis(SpaceNavigatorAxes axisToInvert)
{
	if(axisToInvert < 7 && axisToInvert > 0)
	{
		if(_invertAxes[axisToInvert - 1] == 1.0)
			_invertAxes[axisToInvert - 1] = -1.0;
		else
			_invertAxes[axisToInvert -1] = 1.0;
	}
}

void VRPN_CALLBACK SpaceNavigatorClient::_handleButtons(void *, vrpn_BUTTONCB buttonData)
{
	// read button data
	_spacenavigator->buttons[buttonData.button] = buttonData.state ? true:false;

	if(_spacenavigator->_defaultButtonBehaviour)
	{
		if(_spacenavigator->buttons[0])
			_spacenavigator->switchMode();

		if(_spacenavigator->buttons[1])
			_spacenavigator->switchDomination();
	}
	
	_spacenavigator->_unconsumedData = true;
}

void VRPN_CALLBACK SpaceNavigatorClient::_handleAnalogs(void *, vrpn_ANALOGCB analogData)
{
	// read data from the analog axes
	// range [0, 1]
	// set some values to 0 due to mode

	// normal mode (translation and rotation)
	if(_spacenavigator->_mode == TRANSROT)
	{
		if(_spacenavigator->_dominating)
		{
			double max = analogData.channel[0];
			int index = 0;
			for(int i = 1; i < 6; i++)
			{
				if(abs(analogData.channel[i]) > abs(max))
				{
					index = i;
					max = analogData.channel[i];
				}
			}
			_spacenavigator->x = _spacenavigator->y = _spacenavigator->z = _spacenavigator->rx = _spacenavigator->ry = _spacenavigator->rz = 0;
			switch(index)
			{
			case 0: _spacenavigator->x = max  * _spacenavigator->_invertAxes[0]; break;
			case 2: _spacenavigator->y = max  * _spacenavigator->_invertAxes[1]; break;
			case 1: _spacenavigator->z = max  * _spacenavigator->_invertAxes[2]; break;
			case 3: _spacenavigator->rx = max  * _spacenavigator->_invertAxes[3]; break;
			case 5: _spacenavigator->ry = max  * _spacenavigator->_invertAxes[4]; break;
			case 4: _spacenavigator->rz = max  * _spacenavigator->_invertAxes[5]; break;
			}
		}
		else
		{
			_spacenavigator->x = analogData.channel[0] * _spacenavigator->_invertAxes[0];
			_spacenavigator->y = analogData.channel[2] * _spacenavigator->_invertAxes[1];
			_spacenavigator->z = analogData.channel[1] * _spacenavigator->_invertAxes[2];
			_spacenavigator->rx = analogData.channel[3] * _spacenavigator->_invertAxes[3];
			_spacenavigator->ry = analogData.channel[5] * _spacenavigator->_invertAxes[4];
			_spacenavigator->rz = analogData.channel[4] * _spacenavigator->_invertAxes[5];
		}
	}

	// translation only mode
	else if(_spacenavigator->_mode == TRANS)
	{
		if(_spacenavigator->_dominating)
		{
			double max = analogData.channel[0];
			int index = 0;
			for(int i = 1; i < 3; i++)
			{
				if(abs(analogData.channel[i]) > abs(max))
				{
					index = i;
					max = analogData.channel[i];
				}
			}
			_spacenavigator->x = _spacenavigator->y = _spacenavigator->z = _spacenavigator->rx = _spacenavigator->ry = _spacenavigator->rz = 0;
			switch(index)
			{
			case 0: _spacenavigator->x = max  * _spacenavigator->_invertAxes[0]; break;
			case 2: _spacenavigator->y = max  * _spacenavigator->_invertAxes[1]; break;
			case 1: _spacenavigator->z = max  * _spacenavigator->_invertAxes[2]; break;
			}
		}
		else
		{
			_spacenavigator->x = analogData.channel[0] * _spacenavigator->_invertAxes[0];
			_spacenavigator->y = analogData.channel[2] * _spacenavigator->_invertAxes[1];
			_spacenavigator->z = analogData.channel[1] * _spacenavigator->_invertAxes[2];
			_spacenavigator->rx = _spacenavigator->ry = _spacenavigator->rz = 0;
		}
	}

	// rotation only mode
	else if(_spacenavigator->_mode == ROT)
	{
		if(_spacenavigator->_dominating)
		{
			double max = analogData.channel[3];
			int index = 3;
			for(int i = 4; i < 6; i++)
			{
				if(abs(analogData.channel[i]) > abs(max))
				{
					index = i;
					max = analogData.channel[i];
				}
			}
			_spacenavigator->x = _spacenavigator->y = _spacenavigator->z = _spacenavigator->rx = _spacenavigator->ry = _spacenavigator->rz = 0;
			switch(index)
			{
			case 3: _spacenavigator->rx = max  * _spacenavigator->_invertAxes[3]; break;
			case 5: _spacenavigator->ry = max  * _spacenavigator->_invertAxes[4]; break;
			case 4: _spacenavigator->rz = max  * _spacenavigator->_invertAxes[5]; break;
			}
		}
		else
		{
			_spacenavigator->rx = analogData.channel[3] * _spacenavigator->_invertAxes[3];
			_spacenavigator->ry = analogData.channel[5] * _spacenavigator->_invertAxes[4];
			_spacenavigator->rz = analogData.channel[4] * _spacenavigator->_invertAxes[5];
			_spacenavigator->x = _spacenavigator->y = _spacenavigator->z = 0;
		}
	}

	// swap axes if the z-axis is the up axis
	if(_spacenavigator->_upAxis == Z)
	{
		double temp = _spacenavigator->y;
		_spacenavigator->y = _spacenavigator->z;
		_spacenavigator->z = temp;

		temp = _spacenavigator->ry;
		_spacenavigator->ry = _spacenavigator->rz;
		_spacenavigator->rz = temp;
	}
	
	_spacenavigator->_unconsumedData = true;
	
	#ifdef SPACENAVIGATOR_DEBUG_OUTPUT
	std::cout << "Translation: { " << _spacenavigator->x << ", " << _spacenavigator->y << ", " << _spacenavigator->z << " }" << std::endl;
	std::cout << "Rotation: { " << _spacenavigator->rx << ", " << _spacenavigator->ry << ", " << _spacenavigator->rz << " }" << std::endl;
	#endif // SPACENAVIGATOR_DEBUG_OUTPUT
}

SpaceNavigatorClient* SpaceNavigatorClient::Instance()
{
	if(_spacenavigator == 0)
		_spacenavigator = new SpaceNavigatorClient();
	
	return _spacenavigator;
}
