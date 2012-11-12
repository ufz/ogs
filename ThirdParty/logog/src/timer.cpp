 
/* 
 * \file timer.cpp
 */

#include "logog.hpp"

namespace logog {

	Timer::Timer()
	{
		m_fStartTime = 0.0f;

#ifdef LOGOG_FLAVOR_WINDOWS
		LARGE_INTEGER TicksPerSecond;
		QueryPerformanceFrequency( &TicksPerSecond );
		m_fTicksPerMicrosecond = (DOUBLE)TicksPerSecond.QuadPart * 0.000001;
#endif
		Set( 0.0f );
	}

//! [TimerGet]
	logog::LOGOG_TIME Timer::Get()
	{
#ifdef LOGOG_FLAVOR_WINDOWS
		LARGE_INTEGER liTime;
		QueryPerformanceCounter( &liTime );

		double dusec;
		dusec =( liTime.QuadPart / m_fTicksPerMicrosecond );

		return ( dusec / 1000000.0f ) - m_fStartTime;
#endif

#ifdef LOGOG_FLAVOR_POSIX
#ifdef LOGOG_TARGET_PS3
		LOGOG_PS3_GET_TIME;
#else // LOGOG_TARGET_PS3
		// General Posix implementation
		timeval tv;
		gettimeofday( &tv, 0 );
		return ((double) (tv.tv_sec) + ((double)(tv.tv_usec ) / 1000000.0 ) - m_fStartTime);
#endif // LOGOG_TARGET_PS3
#endif
	}
//! [TimerGet]

	void Timer::Set( LOGOG_TIME time )
	{
		m_fStartTime = time + Get();
	}

	Timer &GetGlobalTimer()
	{
		Statics *pStatic = &Static();

#ifdef LOGOG_INTERNAL_DEBUGGING
		if ( pStatic == NULL )
			LOGOG_INTERNAL_FAILURE;
#endif

		if ( pStatic->s_pTimer == NULL )
			pStatic->s_pTimer = new Timer();

		return *(pStatic->s_pTimer );
	}

	void DestroyGlobalTimer()
	{
		Statics *pStatic = &Static();
		Timer *pGlobalTimer = pStatic->s_pTimer;

		if ( pGlobalTimer != NULL )
			delete pGlobalTimer;

		pStatic->s_pTimer = NULL;
	}

}

