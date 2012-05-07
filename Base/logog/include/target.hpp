/**
 * \file target.hpp Abstractions representing logging outputs.
 */

#ifndef __LOGOG_TARGET_HPP_
#define __LOGOG_TARGET_HPP_

namespace logog
{
/** A target is abstraction representing an output stream from logog.  cerr, cout, and syslog, and other logging formats are supported.
 ** Targets do not validate their received data.  Their only job is to render it on a call to Receive() to their supported target
 ** type.  Targets should generally make sure to handle calls on multiple threads -- they should make sure to avoid overlapping
 ** outputs from multiple threads correctly.  This base class handles this serialization in the Receive() function.
 ** Children of this class are expected to implement the Output() function to do the actual output.
 **/
class Target : public TopicSink
{
    friend class LogBuffer;
public :
    Target();
	virtual ~Target();

    /** Sets the current formatter for this target. */
    void SetFormatter( Formatter &formatter );

    /** Returns a reference to the current formatter for this target. */
    Formatter &GetFormatter() const;

    /** All targets must implement the Output function.  This function outputs the provided string to the
      * output that the target represents.
      * \return Zero if no error has occurred; an error code (generally propagated 
      * from the operating system) if an error has occurred.
      */
    virtual int Output( const LOGOG_STRING &data ) = 0;

    /** Receives a topic on behalf of this target.  A mutex prevents race conditions from occurring when 
     ** multiple threads attempt to write to this target at the same time.
     ** \return Zero if no error has occurred; an error code (generally propagated from the operating
     ** system) if an error has occurred.
     */
    virtual int Receive( const Topic &topic );

	/** Does this target want its formatter to null terminate its strings? */
	bool GetNullTerminatesStrings() const { return m_bNullTerminatesStrings; }

	/** Tells this target whether to request its formatter to null terminate its strings. */
	void SetNullTerminatesStrings(bool val) { m_bNullTerminatesStrings = val; }


protected:
    /** A pointer to the formatter used for this output. */
    Formatter *m_pFormatter;
    /** A mutex on the Receive() function. */
    Mutex m_MutexReceive;
	/** Does this target cause its formatter to null-terminate its output strings?  File and buffer outputs 
	 ** don't require null terminated strings, but line outputs do.
	 **/
	bool m_bNullTerminatesStrings;
};

/** A target representing the cerr stream. */
class Cerr : public Target
{
    virtual int Output( const LOGOG_STRING &data );
};

/** A target representing the cout stream. */
class Cout : public Target
{
    virtual int Output( const LOGOG_STRING &data );
};

/** A target representing the debugger stream on Win32 targets.  This only logs information
  * on Windows like platforms.
  */
class OutputDebug : public Target
{
    virtual int Output( const LOGOG_STRING &data );
};

/** A LogFile renders received messages to a file.  Provide the name of the file to be rendered to as
 * a parameter to the construction of the LogFile() object.  Destroying a LogFile object will cause the
 * output file to be closed.  LogFile objects always append to the output file; they do not delete the previous
 * log file.
 */
class LogFile : public Target
{
public:
    /** Creates a LogFile object.
     * \param sFileName The name of the file to be created.
     * Since file names do not support Unicode on most systems, there is no option to create 
     * a filename with a LOGOG_CHAR. 
     */
    LogFile(const char *sFileName);

    /** Closes the log file. */
    virtual ~LogFile();

    /** Opens the log file on first write. */
    virtual int Open();

	/** This function makes a guess as to the correct BOM for this file, and attempts
	 ** to write it into the file.  It does this by considering the size of LOGOG_CHAR
	 ** as well as considering the current endianness of this system.  This guess
	 ** may be incorrect.
	 **/
	virtual void WriteUnicodeBOM();

    /** Writes the message to the log file. */
    virtual int Output( const LOGOG_STRING &data );

	/** Should a Unicode BOM be written to the beginning of this log file, if the log file
	 * was previously empty?  By default a BOM is written to a log file if LOGOG_UNICODE
	 * is enabled. */
	bool m_bWriteUnicodeBOM;

protected:
    char *m_pFileName;
    bool m_bFirstTime;
	bool m_bOpenFailed;
    FILE *m_pFile;

	/** Does the actual fwrite to the file.  Call Output() instead to handle error conditions better. */
	virtual int InternalOutput( size_t nSize, const LOGOG_CHAR *pData );

private:
    LogFile();
};

/** A buffering target.  Stores up to a fixed buffer size of output and then renders that output to another
  * target.  Can be used for buffering log output in memory and then storing it to a log file upon program completion.
  * To use, create another target (such as a LogFile) and then create a LogBuffer, providing the other target
  * as a parameter to the creation function.
  */
class LogBuffer : public Target
{
public:
    LogBuffer( Target *pTarget = NULL,
               size_t s = LOGOG_DEFAULT_LOG_BUFFER_SIZE );
    virtual ~LogBuffer();

    /** Changes the current rendering target.  NOTE: This function does no locking on either the target or
     * this object.  Program accordingly.
     */
    virtual void SetTarget( Target &t );

    /** Inserts a range of LOGOG_CHAR objects into this buffer.  The characters should consist of a null-terminated
     * string of length size.  Providing anything else as input creates undefined behavior.
     */
    virtual int Insert( const LOGOG_CHAR *pChars, size_t size );

    /** Dumps the current contents of the buffer to the output target. */
    virtual int Dump();

    virtual int Output( const LOGOG_STRING &data );

protected:
    virtual void Allocate( size_t size );

    virtual void Deallocate();

    /** The non-changing pointer to the basic buffer. */
    LOGOG_CHAR *m_pStart;
    /** The current write offset into the buffer. */
    LOGOG_CHAR *m_pCurrent;
    /** The position in the buffer after which no data may be written. */
    LOGOG_CHAR *m_pEnd;
    /** The size of the buffer in LOGOG_CHAR primitives. */
    size_t m_nSize;
    /** A pointer to the target to which the buffer will be rendered upon calling Dump(). */
    Target *m_pOutputTarget;
};

}

#endif // __LOGOG_TARGET_HPP_
