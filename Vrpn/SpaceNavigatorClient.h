// SpaceNavigatorClient.h
//
//
// Author: Lars Bilke

#ifndef SPACENAVIGATORCLIENT_H
#define SPACENAVIGATORCLIENT_H

// comment this out to suppress debug outputs on the console
//#define SPACENAVIGATOR_DEBUG_OUTPUT

// ** INCLUDES **
#include <vrpn_Button.h>
#include <vrpn_Analog.h>

// ** SpaceNavigatorClient - CLASS **
// This implements the singleton design pattern.
// ******************************
class SpaceNavigatorClient
{	
public:
	// ** ENUMS **
	
	/// @brief Navigation modes: translation and rotation or translation
	/// or rotation only.
	enum SpaceNavigatorMode
	{
		TRANSROT = 0,
		TRANS,
		ROT
	};

	/// @brief SpaceNavigator axes.
	enum SpaceNavigatorAxes
	{
		X = 1,
		Y,
		Z,
		rX,
		rY,
		rZ
	};

	// ** Public member fields **

	/// @brief State of the up to eight buttons.
	bool buttons[8];

	// ** Public member functions **

	/// @brief Returns the singleton of the SpaceNavigator.
	static SpaceNavigatorClient* Instance();

	/// @brief Initializes the SpaceNavigator.
	/// Connects with the server and registers the callback handlers.
	/// @param deviceName Example: "SpaceNav@viswork01.intern.ufz.de"
	/// @param axis The up axis.
	/// It is possible to specify the z-axis as the up-axis.
	/// Default the y-axis is the up-axis.
	void init(const char *deviceName, SpaceNavigatorAxes axis = Y);
	
	/// @brief Returns the translation values.
	void getTranslation(double& retx, double& rety, double& retz);
	
	/// @brief Returns the rotation values
	void getRotation(double& retx, double& rety, double& retz);
	
	/// @brief Updates the translation and rotation values.
	/// Must be called once per frame before getTranslation/Rotation.
	void update();
	
	/// @brief Sets the translation scaling factor.
	void setTranslationFactor(double factor);
	
	/// @brief Sets the rotation scaling factor.
	void setRotationFactor(double factor);

	/// @brief Returns true if Z is the up axis.
	bool getZUpAxis();
	
	/// @brief Sets if Z is the up axis.
	void setZUpAxis(bool zUp);

	/// @brief Enables / disables the axis-domination mode.
	/// Only the axis with the highest value is used.
	void setDomination(bool dominating);
	
	/// @brief Switches axis-domination mode.
	void switchDomination();

	/// @brief Switch navigation mode. See SpaceNavigatorMode.
	void setMode(SpaceNavigatorClient::SpaceNavigatorMode mode);
	
	/// @brief Switch through navigation modes.
	void switchMode();

	/// @brief Sets the default button behavior.
	/// on:	left button  --> switch mode
	///		right button --> switch domination
	/// off: no button behavior
	void setDefaultButtonBehaviour(bool enable);

	/// @brief Inverts the specified axis of type SpaceNavigatorAxes
	void invertAxis(SpaceNavigatorAxes axisToInvert);

protected:
	// ** Protected member functions **
	
	/// @brief The constructor is protected because of the singleton
	/// design pattern.
	SpaceNavigatorClient();

	/// @brief Destructor
	virtual ~SpaceNavigatorClient();
	
	/// @brief Returns the elapsed time since the last function call in ms.
	/// This must be implemented in a subclass.
	virtual int getFrameTime() { return 1; }

	/// @brief Does all the event processing.
	/// This function should be run one time per Frame
	/// e.g. in the glut display function.
	// void mainloop();

	/// @brief Actual values of the rotation
	double rx, ry, rz;

	/// @brief Actual values of the translation
	double x, y, z;
	
	/// @brief This one points to the class itself.
	/// You can use only one SpaceNavigatorClient because it´s static.
	/// This is needed for the callback methods which only
	/// can access static members.
	static SpaceNavigatorClient* _spacenavigator;
	
	/// @brief Is set to true if a vrpn callback was called.
	bool _unconsumedData;

private:
	// ** Private member fields **
	vrpn_Button_Remote *_button;
	vrpn_Analog_Remote *_analog;

	// is domination mode active?
	bool _dominating;

	// mode: 0 - translation and rotation
	SpaceNavigatorMode _mode;

	// default button behavior?
	bool _defaultButtonBehaviour;

	// which axis should be inverted?
	double _invertAxes[6];

	// ** Private member functions **

	/// @brief Callbacks which are called whenever a button is pressed
	/// or the SpaceNavigator is moved.
	/// Callbacks as class members have to be static.
	static void VRPN_CALLBACK _handleButtons(void *, vrpn_BUTTONCB buttonData);
	static void VRPN_CALLBACK _handleAnalogs(void *, vrpn_ANALOGCB analogData);
	
	// which is the up-axis (y - normal, z - from the gis world)
	SpaceNavigatorAxes _upAxis;
	
	// The translation factor.
	double _translationFactor;
	
	// The rotation factor.
	double _rotationFactor;
	
	// The translation factor for this frame.
	double _frameTranslationFactor;
	
	// The rotation factor for this frame.
	double _frameRotationFactor;

};
#endif // SPACENAVIGATORCLIENT_H
