 
/* 
 * \file node.cpp
 */

#include "logog.hpp"

namespace logog {

	LockableNodesType & LockableNodesType::operator = (const LockableNodesType &other)
	{
		/* This function is only used at shutdown. */
		LockableNodesType::const_iterator it;

		it = other.begin();
		while ( it != other.end())
		{
			this->insert( *it );
			++it;
		}

		return *this;
	}

	LockableNodesType &GetStaticNodes( void ** pvLocation )
	{
		if ( *pvLocation == NULL )
			*pvLocation = new LockableNodesType();

		return *(( LockableNodesType *)( *pvLocation ));

	}

	LockableNodesType &AllNodes()
	{
	    return GetStaticNodes( &(Static().s_pAllNodes) );
	}

	LockableNodesType &AllSubscriberNodes()
	{
		return GetStaticNodes( &(Static().s_pAllSubscriberNodes ) );
	}

	LockableNodesType &AllFilters()
	{
		return GetStaticNodes( &(Static().s_pAllFilterNodes ) );
	}

	LockableNodesType &AllTargets()
	{
		return GetStaticNodes( &(Static().s_pAllTargets ) );
	}

	Node::Node()
	{
		AllNodes().insert( this );
	}

	Node::~Node()
	{
		Clear();
		AllNodes().erase( this );
	}

	void Node::Initialize()
	{
		if ( CanSubscribe() )
		{
			LockableNodesType *pSubscriberNodes = &AllSubscriberNodes();

			{
				ScopedLock sl( *pSubscriberNodes );
				pSubscriberNodes->insert( this );
			}

			/* This branch is taken iff this node can both subscribe and publish */
			if ( CanPublish() )
			{
				LockableNodesType *pFilterNodes = &AllFilters();
				{
					ScopedLock sl( *pFilterNodes );
					pFilterNodes->insert( this );
				}
			}
		}
	}

	bool Node::CanPublish() const
	{
		return true;
	}

	bool Node::CanSubscribe() const
	{
		return true;
	}

	bool Node::CanSubscribeTo( const Node & )
	{
		return CanSubscribe();
	}

	bool Node::IsTopic() const
	{
		return false;
	}

	bool Node::PublishTo( Node &subscriber )
	{
#ifdef LOGOG_INTERNAL_DEBUGGING
		if ( &subscriber == this )
			LOGOG_INTERNAL_FAILURE;
#endif
		bool bWasInserted;

		{
			ScopedLock sl( m_Subscribers );
			bWasInserted = ( m_Subscribers.insert( &subscriber ) ).second;
		}

		if ( bWasInserted )
			subscriber.SubscribeTo( *this );

		return bWasInserted;
	}

	bool Node::PublishToMultiple( LockableNodesType &nodes )
	{
		LockableNodesType::iterator it;

		bool bWasPublished = false;

		nodes.MutexLock();
		it = nodes.begin();

		while ( it != nodes.end() )
		{
			if ( PublishTo( **it ) == true )
				bWasPublished = true;

			it++;
		}

		nodes.MutexUnlock();

		return bWasPublished;
	}

	bool Node::UnpublishTo( Node &subscriber )
	{
#ifdef LOGOG_INTERNAL_DEBUGGING
		if ( &subscriber == this )
			LOGOG_INTERNAL_FAILURE;
#endif
		bool bWasRemoved = false;

		{
			ScopedLock sl( m_Subscribers );
			NodesType::iterator it;

			if ( ( it = m_Subscribers.find( &subscriber) ) != m_Subscribers.end() )
			{
				bWasRemoved = true;
				m_Subscribers.erase( it );
			}
		}

		if ( bWasRemoved )
			subscriber.UnsubscribeTo( *this );

		return bWasRemoved;
	}

	bool Node::UnpublishToMultiple( LockableNodesType &nodes )
	{
		LockableNodesType::iterator it;

		bool bWasUnpublished = false;

		nodes.MutexLock();
		it = nodes.begin();

		while ( it != nodes.end() )
		{
			if ( UnpublishTo( **it ) == true )
				bWasUnpublished = true;

			it++;
		}

		nodes.MutexUnlock();

		return bWasUnpublished;
	}

	bool Node::SubscribeTo( Node &publisher )
	{
#ifdef LOGOG_INTERNAL_DEBUGGING
		if ( &publisher == this )
			LOGOG_INTERNAL_FAILURE;
#endif
		bool bWasInserted;

		{
			ScopedLock sl( m_Publishers );
			bWasInserted = ( m_Publishers.insert( &publisher ) ).second;
		}

		if ( bWasInserted )
			publisher.PublishTo( *this );

		return bWasInserted;
	}

	bool Node::SubscribeToMultiple( LockableNodesType &nodes )
	{
		LockableNodesType::iterator it;

		bool bWasSubscribed = false;

		nodes.MutexLock();
		it = nodes.begin();

		while ( it != nodes.end() )
		{
			if ( SubscribeTo( **it ) == true )
				bWasSubscribed = true;

			it++;
		}

		nodes.MutexUnlock();

		return bWasSubscribed;
	}

	bool Node::UnsubscribeTo( Node &publisher )
	{
#ifdef LOGOG_INTERNAL_DEBUGGING
		if ( &publisher == this )
			LOGOG_INTERNAL_FAILURE;
#endif
		bool bWasRemoved = false;

		{
			ScopedLock sl( m_Publishers );
			NodesType::iterator it;

			if ( ( it = m_Publishers.find( &publisher ) ) != m_Publishers.end() )
			{
				bWasRemoved = true;
				m_Publishers.erase( it );
			}
		}

		if ( bWasRemoved )
			publisher.UnpublishTo( *this );

		return bWasRemoved;
	}

	bool Node::UnsubscribeToMultiple( LockableNodesType &nodes )
	{
		LockableNodesType::iterator it;

		bool bWasUnsubscribed = false;

		nodes.MutexLock();
		it = nodes.begin();

		while ( it != nodes.end() )
		{
			if ( UnsubscribeTo( **it ) == true )
				bWasUnsubscribed = true;

			it++;
		}

		nodes.MutexUnlock();

		return bWasUnsubscribed;
	}

	void Node::Clear()
	{
		{
			ScopedLock sl( m_Publishers );
			m_Publishers.clear();
		}
		{
			ScopedLock sl( m_Subscribers );
			m_Publishers.clear();
		}
	}


	void DestroyNodesList( void **pvList )
	{
		LockableNodesType **ppNodesList = (LockableNodesType **)pvList;

		if ( *ppNodesList == NULL )
			return;

		(*ppNodesList)->clear();
		delete *ppNodesList;
		*ppNodesList = NULL;
	}

	void DestroyAllNodes()
	{
		Statics *pStatics = &Static();

		LockableNodesType *pAllNodes = ( LockableNodesType *)pStatics->s_pAllNodes;

		if ( pAllNodes == NULL )
			return;

		/** Destroy all the node groups, but don't destroy their contents -- we'll do that as the next step. */
		DestroyNodesList( &(pStatics->s_pAllSubscriberNodes ));
		DestroyNodesList( &(pStatics->s_pAllFilterNodes ));
		DestroyNodesList( &(pStatics->s_pAllTargets ));

		/* We have to copy the AllNodes because destroying each node will remove it from AllNodes.  Fortunately
		 * this only happens at shutdown, so we don't have to worry about efficiency.
		 */
		LockableNodesType nodes = *pAllNodes;

		LockableNodesType::iterator it;

		it = nodes.begin();

		while ( it != nodes.end() )
		{
			delete *it;
			it++;
		}

		nodes.clear();

	#ifdef LOGOG_INTERNAL_DEBUGGING
		if ( pAllNodes->size() != 0 )
			cout << "Not all nodes were deleted at shutdown -- memory leak may have occurred" << endl;
	#endif

		pAllNodes->clear(); // just in case
		delete pAllNodes;
		pStatics->s_pAllNodes = NULL;
	}

}

