/**
 * C++ parser for automata format of APEX
 */

#include "../include/minauto/Util.hpp"
#include "../include/minauto/BLAS.hpp"
#include "../include/minauto/Timer.hpp"
#include "../include/minauto/Automaton.hpp"




std::unordered_map<std::string,int> alphabet_table;
std::vector<std::string> alphabet;



