#include "../include/minauto/Util.hpp"
#include "../include/minauto/RandomValueMap.hpp"

double RandomValueMap::operator()(int x) {
  double result;
  std::unordered_map<int,double>::iterator rand_it(rand_vec.find(x));

  if(rand_it != rand_vec.end()) {
	result = rand_it -> second;
  } else {
	result = getRandomDouble();
	rand_vec.insert(std::pair<int,double>(x,result));
  }
  return result;
}

void RandomValueMap::clear() {
  rand_vec.clear();
}
