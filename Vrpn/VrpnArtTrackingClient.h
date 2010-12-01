#ifndef VRPNARTTRACKINGCLIENT_H
#define VRPNARTTRACKINGCLIENT_H

#include <vrpn_Button.h>
#include <vrpn_Analog.h>
#include <vrpn_Tracker.h>

/// @brief Vrpn client fort ART-Tracking system
class VrpnArtTrackingClient
{
public:
	/// @brief Returns the singleton of this class
	static VrpnArtTrackingClient *Instance();

	/// @brief Initializes and starts the tracking.
	/// @param deviceName The name of the vrpn device @ vrpn server host name
	/// e.g. DTrack@visserv3.intern.ufz.de
	void	StartTracking(const char* deviceName);
	
	/// @brief Stops the tracking.
	void	StopTracking();
	
	/// @brief Returns true if tracking is started
	bool	IsStarted(){return m_bTrackingStarted;}

	/// @brief Calls the vrpn mainloop functions. Must be called once per frame.
	void	MainLoop();

	/// @brief Returns the bodys (head) position.
	void	GetBodyTranslation(double &x, double &y, double &z);
	
	/// @brief Returns the bodys orientation as a quaternion.
	void	GetBodyQuaternion(double &v0, double &v1, double &v2, double &v3);

	/// @brief Returns the flysticks position.
	void	GetFlyTranslation(double &x, double &y, double &z);
	
	/// @brief Returns the flysticks orientation as a quaternion.
	void	GetFlyQuaternion(double &v0, double &v1, double &v2, double &v3);

	/// @brief Returns the analog value of an axis of the flysticks yellow little joystick.
	/// @param index The axis.
	double	GetAnalogData(int index);
	
	/// @brief Returns if the button with the given index of the flystick is pressed.
	bool	GetButtonData(int index);

protected:
	
	/// @brief The constructor is protected because of the singleton
	/// design pattern.
	VrpnArtTrackingClient();
	
	/// @brief Destructor.
	~VrpnArtTrackingClient();
	
	/// @brief This one points to the class itself.
	/// You can use only one VrpnArtTrackingClient because it´s static.
	/// This is needed for the callback methods which only
	/// can access static members.
	static VrpnArtTrackingClient *m_pInstance;

	// Is the tracker initialized ?
	bool	m_bTrackingStarted;

	// Tracking values
	double	m_dBodyQuaternion[4];
	double	m_dBodyTranslation[3];

	// Flystick
	double	m_dFlyQuaternion[4];
	double	m_dFlyTranslation[3];

	// Analogs
	double	m_dAnalogData[10];

	// Buttons
	bool	m_bButtonData[10];

	// VRPN related stuff
	vrpn_Analog_Remote *m_pvrpnAnalog;
	vrpn_Tracker_Remote *m_pvrpnTracker;
	vrpn_Button_Remote *m_pvrpnButtons;

	static void VRPN_CALLBACK CBHandleTracker(void *userdata, const vrpn_TRACKERCB t);
	static void VRPN_CALLBACK CBHandleAnalogs(void *userdata, vrpn_ANALOGCB analogData);
	static void VRPN_CALLBACK CBHandleButtons(void *userdata, vrpn_BUTTONCB buttonData);
};

#endif // VRPNARTTRACKINGCLIENT_H
