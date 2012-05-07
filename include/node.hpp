/**
 * \file node.hpp Base class for higher-level logog objects.
 */

#ifndef __LOGOG_NODE_HPP__
#define __LOGOG_NODE_HPP__

namespace logog
{
/** The base class for most high-level logog objects.  Represents the publisher-subscriber
 ** model within logog. */
class Node;

/** An aggregation of nodes.  Internally, we choose a set representation because we want to be able to traverse
 ** the aggregation quickly while still looking up entries quickly.
 **/
typedef LOGOG_SET< Node *, std::less< Node * >, Allocator< Node * > > NodesType;

/** A type that double inherits from NodesType and Mutex.  A lockable NodesType.  Handles the copy
 ** case correctly.
 **/
class LockableNodesType : public NodesType, public Mutex
{
public:
    /** A LockableNodesType shouldn't copy the internal Mutex when it is copied, but it
     ** should copy all the internal nodes.
     **/
    LockableNodesType & operator = (const LockableNodesType &other);
};

extern LockableNodesType &GetStaticNodes( void ** pvLocation );

/** Returns a reference to the global nodes group.  Allocates a new global node group if one does not already
 ** exist.
 */
extern LockableNodesType &AllNodes();

/** Returns a reference to the group of nodes that are capable of subscribing.  Allocates a new global subscriber
 ** node group if one does not already exist.
 */
extern LockableNodesType &AllSubscriberNodes();

/** Returns a reference to the group of nodes that are capable of both subscribing as well as publishing.  Allocates a new global subscriber
 ** node group if one does not already exist.
 */
extern LockableNodesType &AllFilters();

/** Returns a reference to the group of nodes that represent terminals in the graph, i.e. nodes that can't publish. */
extern LockableNodesType &AllTargets();

class Node : public Object
{
public:

    /** All nodes self-register as part of the all-nodes database. */
    Node();

    ~Node();

    /** Call this function immediately after creating a node (or any of the children of the node class.)  This function currently
     ** registers the node as part of the list of subscriber nodes, if this node may in fact subscribe.
     ** If this node is capable of subscribing at all, then this function registers this node as a possible subscriber.
     ** Doing this helps to keep down the number of nodes we search, when we are determining which nodes a new node
     ** might subscribe to.  We have to do this registration as a second step, after the node is completely
     ** initialized, as subscriberness is determined late in initialization.
     **/
    virtual void Initialize();

    /** Can a node send notifications?  By default they can; later subclasses may not be able to. */
    virtual bool CanPublish() const;
    /** Can a node receive notifications?  By default they can; later subclasses may not be able to. */
    virtual bool CanSubscribe() const;
    /** Is this node interested in receiving notifications from another topic? */
    virtual bool CanSubscribeTo( const Node & );

    /** In order to avoid bringing in a bunch of RTTI stuff, we permit nodes to be asked whether they're topics or not */
    virtual bool IsTopic() const;

    /** Causes this node to begin publishing events to the given subscriber.
     ** \param subscriber The node to receive published events
     ** \return true if the request was successful, false if the subscriber was already subscribed
     **/
    virtual bool PublishTo( Node &subscriber );

    /** Causes this node to attempt to publish to some other nodes. */
    virtual bool PublishToMultiple( LockableNodesType &nodes );

    /** Causes this node to stop publishing events to this subscriber.
     ** \param subscriber The node to stop receiving events
     ** \return true if successful, false if the subscriber was not being published to in the first place
     **/
    virtual bool UnpublishTo( Node &subscriber );

    /** Causes this node to attempt to unpublish to some other nodes. */
    virtual bool UnpublishToMultiple( LockableNodesType &nodes );

    /** Causes this node to start receiving events from the given publisher.
     ** \param publisher The node to start receiving events from
     ** \return true if successful, false if the publisher was already subscribed
     **/
    virtual bool SubscribeTo( Node &publisher );

    /** Causes this node to attempt to subscribe to some other nodes. */
    virtual bool SubscribeToMultiple( LockableNodesType &nodes );


    /** Causes this node to unsubscribe from the given publisher's events.
    ** \param publisher The publisher to unsubscribe from
    ** \return true if successful, false if the node was already unsubscribed
    **/
    virtual bool UnsubscribeTo( Node &publisher );

    /** Causes this node to attempt to unsubscribe to some other nodes. */
    virtual bool UnsubscribeToMultiple( LockableNodesType &nodes );

    void Clear();


    /** A pointer to any custom data you need to store for a node. */
    void *m_pUserData1;

    /** A pointer to any custom data you need to store for a node. */
    void *m_pUserData2;

protected:
    /** A bunch of nodes that are interested in what this node has to report. */
    LockableNodesType	m_Subscribers;

    /** A bunch of nodes that this node interested in hearing from. */
    LockableNodesType	m_Publishers;
};

extern void DestroyNodesList( void **pvList );

/** Destroys all nodes currently recorded.  This happens at shutdown time.  NOTE!  If you have allocated
 ** your own logog items and free them yourself AFTER this call, exciting crashes will occur.
 **/
extern void DestroyAllNodes();

}

#endif // __LOGOG_NODE_HPP_
