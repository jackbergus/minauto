#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <sys/time.h>
using namespace std;

#include "../include/minauto/Timer.hpp"

void Error(const char* msg) {
  printf("%s\n",msg);
  abort();
}


/*********************************************************************/
//constructors and destructors
/*********************************************************************/
Timer::Timer(const Timer &rhs) { *this = rhs; }

/*********************************************************************/
//operators
/*********************************************************************/
const Timer &Timer::operator = (const Timer &rhs)
{
  start = rhs.start;
  time = rhs.time;
  return *this;
}

/*********************************************************************/
//return true if the timer is running and false otherwise
/*********************************************************************/
bool Timer::Running() const { return (start != -1); }

/*********************************************************************/
//read the current time. the timer must be stopped.
/*********************************************************************/
double Timer::Read() const
{
  if(start != -1) {
    Error("ERROR: trying to read a running timer ...\n");
  }
  return time;
}

/*********************************************************************/
//forward a timer by the supplied amount of time
/*********************************************************************/
void Timer::Forward(const Timer &arg)
{
  if(arg.start != -1) {
    Error("ERROR: trying to forward by a running timer ...\n");
  }
  if(start != -1) {
    Error("ERROR: trying to forward a running timer ...\n");
  }
  time += arg.time;
}

/*********************************************************************/
//start the timer
/*********************************************************************/
void Timer::Start()
{
  if(start != -1) {
    Error("ERROR: trying to start a running timer ...\n");
  }
#ifdef WIN32
  start = static_cast<double>(::time(NULL));
#else
  struct timeval tv;
  if(gettimeofday(&tv,NULL) == -1) {
    Error("ERROR: could not get current time for starting timer ...\n");
  }
  start = tv.tv_sec + (tv.tv_usec / 1000000.0);
#endif //WIN32
}

/*********************************************************************/
//stop the timer
/*********************************************************************/
void Timer::Stop()
{  
  if(start == -1) {
	Error("WARNING: Timer already stopped\n");
	return;  
  }
#ifdef WIN32
  double currTime = static_cast<double>(::time(NULL));
#else
  struct timeval tv;
  if(gettimeofday(&tv,NULL) == -1) {
    Error("ERROR: could not get current time for stopping timer ...\n");
  }
  double currTime = tv.tv_sec + (tv.tv_usec / 1000000.0);
#endif //WIN32
  time += (currTime - start);
  start = -1;
}


/*********************************************************************/
//end of Timer.cpp
/*********************************************************************/
