 /* 
 * \file statics.cpp
 */

#include "logog.hpp"

namespace logog {

	Statics::Statics()
	{
		s_pAllNodes = NULL;
		s_pAllSubscriberNodes = NULL;
		s_pAllFilterNodes = NULL;
		s_pAllTargets = NULL;
		s_pTimer = NULL;
		s_pDefaultFormatter = NULL;
		s_pDefaultFilter = NULL;
		s_pStringSearchMutex = NULL;
		s_pMessageCreationMutex = NULL;
		s_pfMalloc = NULL;
		s_pfFree = NULL;
		s_pSelf = this;
		s_nSockets = 0;
	}

	void Statics::Reset()
	{
		DestroyGlobalTimer();
		DestroyDefaultFormatter();
		s_pDefaultFilter = NULL; // This will be destroyed on the next step
		DestroyAllNodes();
		DestroyStringSearchMutex();
		DestroyMessageCreationMutex();
		s_pfMalloc = NULL;
		s_pfFree = NULL;
		s_nSockets = 0;
	}

	Statics::~Statics()
	{
		Reset();
	}

	Statics s_Statics;

	Statics &Static()
	{
		return s_Statics;
	}

	void DestroyStatic()
	{
		s_Statics.~Statics();
	}

}

