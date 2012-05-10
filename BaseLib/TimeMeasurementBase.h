#ifndef TIMEMEASUREMENT_H
#define TIMEMEASUREMENT_H

namespace BaseLib {

class TimeMeasurementBase
{
public:
	virtual void start () = 0;
	virtual void stop () = 0;
	virtual double elapsed () = 0;
	virtual ~TimeMeasurementBase () {};
};

} // end namespace BaseLib

#endif
