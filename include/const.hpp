/**
 * \file const.hpp Constants.
 */

#ifndef __LOGOG_CONST_HPP__
#define __LOGOG_CONST_HPP__

#ifndef LOGOG_FORMATTER_MAX_LENGTH
/** The maximum length of a single line that a formatter may output, in LOGOG_CHAR units. */
#define LOGOG_FORMATTER_MAX_LENGTH ( 1024 * 16 )
#endif

#ifndef LOGOG_DEFAULT_LOG_BUFFER_SIZE
/** The default size of a RingBuffer object for buffering outputs. */
#define LOGOG_DEFAULT_LOG_BUFFER_SIZE ( 4 * 1024 * 1024 )
#endif


/** \addtogroup levelsettings Level Settings
 ** These are level settings for logog.  These settings are valid for the LOGOG_LEVEL compilation flag.  To enable
 ** all logog messages, use the compilation flag -DLOGOG_LEVEL=LOGOG_LEVEL_ALL.  To disable all logog messages
 ** (and effectively remove logog code from your executable) use the compilation flag -DLOGOG_LEVEL=LOGOG_LEVEL_NONE.
 ** \sa LOGOG_LEVEL
 ** @{
 **/
//! [Level Constants]
#define LOGOG_LEVEL_NONE		0
#define LOGOG_LEVEL_EMERGENCY	8
#define LOGOG_LEVEL_ALERT		16
#define LOGOG_LEVEL_CRITICAL	24
#define LOGOG_LEVEL_ERROR		32
#define LOGOG_LEVEL_WARN		40
#define LOGOG_LEVEL_WARN1		48
#define LOGOG_LEVEL_WARN2		56
#define LOGOG_LEVEL_WARN3		64
#define LOGOG_LEVEL_INFO		72
#define LOGOG_LEVEL_DEBUG		80
#define LOGOG_LEVEL_ALL			88
//! [Level Constants]

#define LOGOG_LEVEL_TYPE		int

#ifndef LOGOG_LEVEL
#define LOGOG_LEVEL LOGOG_LEVEL_DEBUG
#endif

/** @} */

/** \addtogroup topicbitstype Topic Bits Type
  * Bit flags representing whether a topic cares about a specific field or not.  1 = care, 0 = don't care.
  * @{
  */
//! [Topic Bits]
typedef enum
{
	TOPIC_LEVEL_FLAG =			0x01,
	TOPIC_LINE_NUMBER_FLAG =	0x02,
	TOPIC_FILE_NAME_FLAG =		0x04,
	TOPIC_GROUP_FLAG =			0x08,
	TOPIC_CATEGORY_FLAG =		0x10,
	TOPIC_MESSAGE_FLAG =		0x20,
	TOPIC_TIMESTAMP_FLAG =		0x40,
	/** Bits 0 through TOPIC_COUNT turned on */
	TOPIC_ALL = 0x7f
} TopicBitsType;
//! [Topic Bits]

/** @} */

typedef int TOPIC_FLAGS;

//! [Topic Offsets]
/** Offsets within the m_vIntProps and m_vStringProps arrays for this topic. */
typedef enum
{
	TOPIC_LEVEL = 0,
	TOPIC_LINE_NUMBER = 1,
	/** This must be the number of integer fields. */
	TOPIC_INT_COUNT = 2,

	TOPIC_FILE_NAME = 0,
	TOPIC_GROUP = 1,
	TOPIC_CATEGORY = 2,
	TOPIC_MESSAGE = 3,
	/** This must be the number of string fields for this topic. */
	TOPIC_STRING_COUNT = 4
} TopicOffsetType;
//! [Topic Offsets]

#endif // __LOGOG_CONST_HPP__
