/**
 * \file timer.hpp Time management.
 */

#ifndef __LOGOG_TIMER_HPP__
#define __LOGOG_TIMER_HPP__

namespace logog
{
/** A value for a high resolution timer on this platform.  Time representations are in seconds. */
typedef double LOGOG_TIME;

/** A high-resolution timer.  Reports in seconds. */
class Timer : public Object
{
public:
    Timer();

    /** Returns the offset from the time since the creation of the timer, or the time set by the most
     ** recent Set() call.  Time is assumed to be a value in LOGOG_TIME seconds.
     ** \sa LOGOG_TIME
     **/
    LOGOG_TIME Get();

    /** Sets the current time for this timer. */
    void Set( LOGOG_TIME time );

protected:
#ifdef LOGOG_FLAVOR_WINDOWS
    /** Windows only.  Stores the number of high resolution timer ticks per second. */
    double m_fTicksPerMicrosecond;
#endif
    /** Zero, if no calls to Set() have been made; else the value of the previous call to Set(). */
    LOGOG_TIME m_fStartTime;
};

extern Timer &GetGlobalTimer();
extern void DestroyGlobalTimer();

}

#endif // __LOGOG_TIMER_HPP_
