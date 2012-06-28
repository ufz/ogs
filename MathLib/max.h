/**
 * \file max.h
 *
 * Created on 2010-04-20 by Thomas Fischer
 */

#ifndef MAX_H_
#define MAX_H_

/**
 * max returns the maximum of its arguments
 */
template<class T> T max(const T& arg0, const T& arg1)
{
  if (arg0 < arg1) return arg1;
  return arg0;
}

#endif /* MAX_H_ */
