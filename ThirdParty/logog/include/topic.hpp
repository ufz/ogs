/**
 * \file topic.hpp Topics -- subjects of interests for publishers and
 * subscribers to communicate about
 */

#ifndef __LOGOG_TOPIC_HPP__
#define __LOGOG_TOPIC_HPP__

namespace logog
{

/** A subject that nodes can choose to discuss with one another.
 ** Subscribers generally have very general topics, while publishers generally have very specific topics.
 **/
class Topic : public Node
{
    friend class TopicLevel;
    friend class TopicGroup;

public:
    /** Creates a topic.  Note the defaults for creating a topic -- these defaults are equivalent to "no setting"
     ** for those fields.
     **/
    Topic( const LOGOG_LEVEL_TYPE level = LOGOG_LEVEL_ALL,
           const LOGOG_CHAR *sFileName = NULL,
           const int nLineNumber = 0,
           const LOGOG_CHAR *sGroup = NULL,
           const LOGOG_CHAR *sCategory = NULL,
           const LOGOG_CHAR *sMessage = NULL,
           const double dTimestamp = 0.0f );

    /** Topics are always topics.  We use this to avoid any RTTI dependence. */
    virtual bool IsTopic() const;

    /** Causes this topic to publish another topic to all its subscribers.
     ** \return 0 if successful, non-zero if this topic failed to send the publication to all subscribers */
    virtual int Send( const Topic &node );

    /** Causes this topic to publish itself to all its subscribers. */
    virtual int Transmit();

    /** Permits this node to receive a publication from another node, and act upon it.
     ** \param node The node constituting the publication
     ** \return 0 if successful, non-zero if this node failed to process the publication
     **/
    virtual int Receive( const Topic &node );

    /** Is this topic interested in receiving notifications from another topic?  This function implements
     ** a generic (slow) test that should work for all topic types.  This function only checks fields
     ** that have previously been set on this topic -- fields that have not been set will not limit this
     ** topic's ability to subscribe.  If any of the previously set fields does not "match" the other topic,
     ** this function will return false.  The matching function behaves slightly differently from field to
     ** field.
     ** - In the topic level case, this function rejects a publisher with a lower topic level than our
     ** own.
     ** - In the group, category, file name and message case, this function rejects a publisher if our
     ** own group/category/file name or message cannot be found as a substring in the possible publisher.
     ** This functionality permits very simple pattern matching functionality (i.e. show me all the message
     ** lines that have the word "upload" in them, regardless of their log level.)
     ** - In the line number case, this function rejects a publisher unless the line number matches exactly.
     ** - In the timestamp case, this function rejects a publisher if their timestamp is before our own.
     ** \param otherNode The topic which we are considering subscribing to
     **/
    virtual bool CanSubscribeTo( const Node &otherNode );

    virtual bool CanSubscribeCheckTopic( const Topic &other );

    /** Causes this topic to begin publishing events to the given subscriber.
     ** \param subscriber The node to receive published events
     ** \return true if the request was successful, false if the subscriber was already subscribed
     **/
    virtual bool PublishTo( Node &subscriber );



    /** Formats the message in this topic given a sprintf-style set of arguments.
     ** This function can be used to set the current message in this topic to a string with a variable number of parameters.
     **/
    virtual void Format( const LOGOG_CHAR *cFormatMessage, ... );

    const LOGOG_STRING &FileName() const;
    void FileName( const LOGOG_STRING &s );

    const LOGOG_STRING &Message() const;
    void Message( const LOGOG_STRING &s );

    const LOGOG_STRING &Category() const;
    void Category( const LOGOG_STRING &s );

    const LOGOG_STRING &Group() const;
    void Group( const LOGOG_STRING &s );

    int LineNumber() const;
    void LineNumber( const int num );

    LOGOG_LEVEL_TYPE Level() const;
    void Level( LOGOG_LEVEL_TYPE level );

    LOGOG_TIME Timestamp() const;
    void Timestamp( const LOGOG_TIME t );

    TOPIC_FLAGS GetTopicFlags() const;

protected:
    /** An array (not an STL vector) of string properties for this topic. */
    LOGOG_STRING m_vStringProps[ TOPIC_STRING_COUNT ];
    /** An array (not an STL vector) of integer properties for this topic. */
    int m_vIntProps[ TOPIC_INT_COUNT ];
    /** The time associated with this topic.  Usually this field is updated when a topic is triggered.  Times need not be associated
     ** with a particular topic, in which case this value is zero.
     ** */
    LOGOG_TIME m_tTime;
    /** A bitfield representing the "important" fields in this topic.  Not all fields are considered to contain important information
     ** all the time.  A logical OR of the TOPIC_*_FLAG fields.
     ** \sa TopicBitsType
     **/
    TOPIC_FLAGS m_TopicFlags;
};

/** A topic that permits both publishing as well as subscribing.  This class is functionally same as a Topic; we've added it
 ** as a class for clarity when referring to different topic types.  Filters should be instantiated only after outputs are
 ** instantiated, as they automatically search for and publish to targets.  If you instantiate a filter before you
 ** instantiate a target, you will need to call PublishTo( theTarget ) yourself before using the target.
 **/
class Filter : public Topic
{
public:
    Filter( const LOGOG_LEVEL_TYPE level = LOGOG_LEVEL_ALL,
            const LOGOG_CHAR *sFileName = NULL,
            const int nLineNumber = 0,
            const LOGOG_CHAR *sGroup = NULL,
            const LOGOG_CHAR *sCategory = NULL,
            const LOGOG_CHAR *sMessage = NULL,
            const double dTimestamp = 0.0f );
};

/** Returns a reference to the unique default filter instantiated with logog. */
extern Filter &GetDefaultFilter();

/** Sets the current reporting level for the default filter.  All messages
  * connected to the filter after this point should obey this default
  * level setting.
  */
void SetDefaultLevel( LOGOG_LEVEL_TYPE level );

/** A topic with the group name being the only field of significance. */
class TopicGroup : public Topic
{
public:
    TopicGroup( const LOGOG_CHAR *sGroup = NULL );

    virtual bool CanSubscribeCheckTopic( const Topic &other );
};

/** A topic with the level being the only field of significance. */
class TopicLevel : public Topic
{
public:
    TopicLevel( const LOGOG_LEVEL_TYPE level );
    virtual bool CanSubscribeCheckTopic( const Topic &other );
};

/** A topic that is also a source. */
class TopicSource : public Topic
{
public:
    TopicSource( const LOGOG_LEVEL_TYPE level = LOGOG_LEVEL_ALL,
                 const LOGOG_CHAR *sFileName = NULL,
                 const int nLineNumber = 0,
                 const LOGOG_CHAR *sGroup = NULL,
                 const LOGOG_CHAR *sCategory = NULL,
                 const LOGOG_CHAR *sMessage = NULL,
                 const double dTimestamp = 0.0f );

    /** Returns false.  Sources do not subscribe. */
    virtual bool SubscribeTo( Node & );

    /** Returns false.  Sources do not unsubscribe. */
    virtual bool UnsubscribeTo( Node & );
    virtual bool CanSubscribe() const;
};

/** A topic that is also a sink. */
class TopicSink : public Topic
{
public:
    virtual bool IsTopic() const;

    /** Sinks do not add themselves to the list of interested subscribers.  That's up to intermediate topics to decide. */
    virtual void Initialize();

    /** Returns false.  Sinks do not publish. */
    virtual bool PublishTo( Node & );

    /** Returns false.  Sinks do not unpublish. */
    virtual bool UnpublishTo( Node & );

    /** Returns false.  Sinks do not publish. */
    virtual bool CanPublish() const;
};
}

#endif // __LOGOG_TOPIC_HPP_
