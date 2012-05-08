/* 
 * \file message.cpp
 */

#include "logog.hpp"

namespace logog
{
	Message::Message( const LOGOG_LEVEL_TYPE level,
             const LOGOG_CHAR *sFileName ,
             const int nLineNumber,
             const LOGOG_CHAR *sGroup,
             const LOGOG_CHAR *sCategory,
             const LOGOG_CHAR *sMessage,
             const double dTimestamp,
			 bool *pbIsCreated ) :
        Checkpoint( level, sFileName, nLineNumber, sGroup, sCategory, sMessage, dTimestamp )
    {
		m_pbIsCreated = pbIsCreated;

		if ( pbIsCreated != NULL )
			*pbIsCreated = true;

        /* Messages are always sources, so there's no need to call Initialize() here */
        // Initialize();

        /* NOTE!  The message is typically assigned to a checkpoint AFTER it's been published.
         * Ergo a message needs to be unpublished from, and published to, all filters
         * iff one of those filters is searching for a substring of that message.
         */
        PublishToMultiple( AllFilters() );
    }

	Message::~Message()
	{
		if ( m_pbIsCreated != NULL )
			*m_pbIsCreated = false;
	}


	bool Message::Republish()
	{
		UnpublishToMultiple( AllFilters() );
		return PublishToMultiple( AllFilters() );
	}

	Mutex &GetMessageCreationMutex()
	{
		Statics *pStatic = &Static();

		if ( pStatic->s_pMessageCreationMutex == NULL )
			pStatic->s_pMessageCreationMutex = new Mutex();

		return *( pStatic->s_pMessageCreationMutex );
	}

	void DestroyMessageCreationMutex()
	{
		Statics *pStatic = &Static();

		if ( pStatic->s_pMessageCreationMutex != NULL )
		{
			delete pStatic->s_pMessageCreationMutex;
			pStatic->s_pMessageCreationMutex = NULL;
		}
	}
}


