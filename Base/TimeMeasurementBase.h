#ifndef TIMEMEASUREMENT_H
#define TIMEMEASUREMENT_H

class TimeMeasurementBase 
{
public:
	virtual void start () = 0;
	virtual void stop () = 0;
	virtual double elapsed () = 0;
	virtual ~TimeMeasurementBase () {};
};

#endif

