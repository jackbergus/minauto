#ifndef __TIMER_H__
#define __TIMER_H__

/*! \class Timer The class allows measurment of time 
    \ingroup util
*/
class Timer
{
 private:
  double start;
  double time;

 public:
  Timer() { start = -1; time = 0; }
  Timer(const Timer &rhs);  
  const Timer &operator = (const Timer &rhs);

  bool Running() const;
  void Start();
  void Stop();
  void Forward(const Timer &arg);
  double Read() const;
};

#endif //__TIMER_H__
