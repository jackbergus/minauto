#ifndef WORDTREE
#define WORDTREE

#include "Automaton.hpp"

extern std::vector<std::string> alphabet;

struct WordTree {
  // a node consists of 
  // * (1st component) a letter 
  // * (2nd component) a parent-pointer index
  typedef std::pair<int,int> Node;
  std::vector< Node > nodes;
  
  int addNode(int a, int b);

  void getWord(int b, Word& w, bool forward);
};

#endif
