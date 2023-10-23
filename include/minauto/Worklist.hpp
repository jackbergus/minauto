#ifndef WORKLIST
#define WORKLIST

#include <list>

/*
  - a pair of distributions and a word
 */
template<class Vector>
struct DWPair {
  Vector ua;
  Vector ub;
  int node;

  DWPair (const Vector& __ua, const Vector& __ub, int __node) : ua(__ua), ub
(__ub), node(__node) {}
  DWPair () {}
};

/*
  - a distributions and a word
 */
template<class Vector>
struct DW {
  Vector u;
  int node;
  int depth;

  DW (const Vector& __u, int __node, int __depth) : u(__u), node(__node), depth(__depth) {}
  DW () {}
};



template<class Entry>
class Worklist {
public:
  Worklist ()  {}

  void add(const Entry& v) {
    worklist.push_back(v);
  }

  void clear() {
    worklist.clear();
  }

  Entry& top() {
    return worklist.front();
  }

  void pop() {
    worklist.pop_front();
  }

  unsigned size () const { return worklist.size(); }
protected:
  std::list<Entry> worklist;
};


#endif
