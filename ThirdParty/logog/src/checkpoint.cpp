 
/* 
 * \file checkpoint.cpp
 */

#include "logog.hpp"

namespace logog {

	Checkpoint::Checkpoint( const LOGOG_LEVEL_TYPE level,
		const LOGOG_CHAR *sFileName ,
		const int nLineNumber,
		const LOGOG_CHAR *sGroup,
		const LOGOG_CHAR *sCategory,
		const LOGOG_CHAR *sMessage,
		const double dTimestamp ) :
	TopicSource( level, sFileName, nLineNumber, sGroup, sCategory, sMessage, dTimestamp )
	{
	}

	int Checkpoint::Send( const Topic &node )
	{
		/* Optionally update our own timestamp before we send on our information */
		if (( m_TopicFlags & TOPIC_TIMESTAMP_FLAG ) != 0 )
			m_tTime = GetGlobalTimer().Get();

		return TopicSource::Send( node );
	}
}

