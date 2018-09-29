 /* 
 * \file api.cpp
 */

#include "logog.hpp"

#include <cstdlib>

namespace logog {

static Mutex s_mInitialization;
static int s_nInitializations = 0;

int Initialize( INIT_PARAMS *params )
{
	s_mInitialization.MutexLock();

	if ( s_nInitializations++ == 0 )
	{
		if ( params == NULL )
		{
			Static().s_pfMalloc = malloc;
			Static().s_pfFree = free;
		}
		else
		{
			if ( params->m_pfMalloc != NULL )
			{
				Static().s_pfMalloc = params->m_pfMalloc;
				Static().s_pfFree = params->m_pfFree;
			}
			else
			{
				Static().s_pfMalloc = malloc;
				Static().s_pfFree = free;
			}
		}

		// Let's allocate a default filter here.
		GetDefaultFilter();

		// Socket::Initialize();
	}

	s_mInitialization.MutexUnlock();

    return 0;
}

int Shutdown( )
{
	s_mInitialization.MutexLock();

	if ( --s_nInitializations == 0 )
	{
		// Socket::Shutdown();

#ifdef LOGOG_DESTROY_STATIC_AREA
		delete &( Static() );
#else
		Static().Reset();
#endif

#ifdef LOGOG_LEAK_DETECTION
		ReportMemoryAllocations();
#endif
	}

	s_mInitialization.MutexUnlock();

    return 0;
}
}

