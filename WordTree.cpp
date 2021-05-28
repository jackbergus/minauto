#include <string>
#include <vector>
#include <list>
#include <iostream>

typedef std::list<int> Word;

extern std::vector<std::string> alphabet;

#include "WordTree.hpp"

int WordTree::addNode(int a, int b) {
	nodes.push_back(Node(a,b));
	return nodes.size() - 1; // index to the added entry
}

void WordTree::getWord(int b, Word& word, bool forward) {
	if(b > nodes.size() && b < 0)
		return;

	if(!forward) {
		for(int current = b; current > 0; current = nodes[current].second) {
			word.push_back(nodes[current].first);
		} 
	} else {
		for(int current = b; current > 0; current = nodes[current].second) {
			word.push_front(nodes[current].first);
		} 
	}
}


