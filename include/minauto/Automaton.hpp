#ifndef AUTOMATON
#define AUTOMATON

#include <list>
#include <string>
#include <unordered_map>
#include <vector>
#include <fstream>
#include "Util.hpp"
#include <BLAS.hpp>

#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>

/* global alphabet table for both automata */
extern std::unordered_map<std::string,int> alphabet_table;
extern std::vector<std::string> alphabet;

typedef std::list<int> Word;

template<typename X>
class Automaton {
public:
  typedef X Weight;
  typedef typename boost::numeric::ublas::vector<X> Vector;
  typedef typename boost::numeric::ublas::zero_vector<X> ZeroVector;

  typedef typename Vector::iterator VectorIterator;
  typedef typename boost::numeric::ublas::vector_range< Vector > VectorRange;
  typedef boost::numeric::interval<X> Interval;
  typedef boost::numeric::ublas::vector<Interval, std::vector<Interval> > IntervalVector;
  typedef boost::numeric::ublas::compressed_matrix<X,boost::numeric::ublas::row_major> Matrix;

  typedef typename boost::numeric::ublas::identity_matrix<X> IdentityMatrix;
  typedef typename boost::numeric::ublas::matrix_range< Matrix > MatrixRange;


  typedef std::unordered_map<int,Matrix> TransitionMatrix;

  typedef std::unordered_map<int,std::vector<unsigned> > TransitionNZ;


  Automaton (int n = 0) : initial(n,1), accepting(n,1), nr_of_states(n), nr_of_transitions(0), accepting_state (0), buffer_length(20000) { initial.clear(); accepting.clear(); }

  // a - b
  Automaton(const Automaton& a, const Automaton& b);
   
  template<typename Y>
  Automaton(const Automaton<Y>& a);

  ~Automaton () {
    alphabet.clear();
    trans.clear();
  }

  void setNrOfStates(unsigned int n) { nr_of_states = n; incoming.resize(n); enabled.resize(n); accepting.resize(n); initial.resize(n);}

  unsigned getNrOfStates() const { return nr_of_states; }
  unsigned getNrOfTransitions() const { return nr_of_transitions; }


  unsigned getNrOfNonzeros();

  void setInitial(int state, X weight = 1.0) { initial(state) = weight; }
  void setAccepting(int state, X weight = 1.0) { accepting(state) = weight; }

  void addTransition(unsigned int source, unsigned int target, std::string sigma, X prob);
  void addTransition(unsigned int source, unsigned int target, int sigma, X prob);
  void addDistribution(int state, int value, X probability);

  void computeMatrices();

  long getMem() const;

  std::vector<std::unordered_set<int> > enabled;
  std::vector<std::unordered_set<int> > incoming;

  TransitionMatrix trans_matrix;
  TransitionNZ trans_nz;

  std::unordered_set<int> accepting_actions;
 
  Vector initial;
  Vector accepting;

  const Vector& getStart(bool forward) const { return forward ? initial : accepting; }
  const Vector& getFinish(bool forward) const { return forward ? accepting : initial; }



  const Vector& getInitial() const { return initial; }
  const Vector& getAccepting() const { return accepting; }



  // compute probability of a given word
  X operator()(const Word& w) const;

  // input in  format
  int parse(const char* filename) {
	  int result = read(filename);
	  if(result) {
		  // this might be an apex file
		  result = readApex(filename);
	  }
	  return result;
  }


  // input in generic format
  int read(const char* filename);

  // output in generic format
  int write(const char* filename);

  // Matlab matrix output
  int writeMatlab(const char* filename);

  // read automaton from APEX file
  int readApex(const char* filename);

  // output into APEX format
  int writeApex(const char* filename);

  // graph output into aiSee format
  // http://www.aisee.com/
  void  aiSee(std::ostream& stream) const;


  void getLetters(const Vector& v, std::unordered_set<int>& result, bool forward) {
	  if(forward)
		  get_enabled(v,result);
	  else
		  get_incoming(v,result);
  }

  void get_enabled(const Vector& v, std::unordered_set<int>& result) const {
  	typedef typename Vector::const_iterator it_t;
	for(it_t it=v.begin(); it!=v.end();++it) {
		if(*it == 0.0) continue;
		unsigned index = it.index();
		const std::unordered_set<int>& s(enabled[index]);
	    result.insert(s.begin(),s.end());
	}
  }

  void get_incoming(const Vector& v, std::unordered_set<int>& result) const {
  	typedef typename Vector::const_iterator it_t;
	for(it_t it=v.begin(); it!=v.end();++it) {
	    if(*it == 0.0) continue;
	    unsigned index = it.index();
	    const std::unordered_set<int>& s(incoming[index]);
	    result.insert(s.begin(),s.end());
	}
   }

  template <class V>
  void operator() (int sigma, V& next, const V& current, bool forward = true) const {
	  typename TransitionMatrix::const_iterator mit = trans_matrix.find(sigma);
	  typename TransitionNZ::const_iterator nit = trans_nz.find(sigma);

	  if( mit == trans_matrix.end() || nit == trans_nz.end()) {
		  next = boost::numeric::ublas::zero_vector<X>(current.size());
	  } else {
		  if(forward) {
			  row_prod(next,current,mit->second,nit->second);
		  } else {
			  col_prod(next,current,mit->second,nit->second);

		  }
	  }
  }

  void reachabilityAnalysis(bool forward) {
	  Vector start ( forward ? initial : accepting);

	  unsigned counter = 0;


	  std::cout << "reachability " << counter << " states " << std::endl;
  }


  class Transition {
  public:
    unsigned int source, target;
    X prob;
    int sigma;
  };

  struct TransitionComp {
    bool operator() (const Transition& t1, const Transition& t2) {
  	return t1.sigma < t2.sigma || ( t1.sigma == t2.sigma && t1.source < t2.source ) ;
    }
  };

protected:      
  unsigned int nr_of_states;
  unsigned int nr_of_transitions;

  std::vector<Transition> trans;
  unsigned accepting_state;
  unsigned buffer_length;


  int readVector(char* input, Vector& vec) const;
  int parseDistribution(char* input, int accepting_state);
  int parseTransition(char* input);
  bool checkRangeOfStateIndex(int) const;
  int writeTransitions(std::ofstream& f);
  int writeVector(const Vector& v, std::ofstream& f);
  int writeMatlabMatrix(std::ofstream& f) ;
  int writeMatlabVector(const Vector& v, std::ofstream& f) ;
};


template<typename X>
template<typename Y>
Automaton<X>::Automaton(const Automaton<Y>& a) :
	enabled(a.enabled),
        incoming(a.incoming),
	initial(a.initial),
	accepting(a.accepting),
	nr_of_states (a.getNrOfStates()),
	nr_of_transitions (a.getNrOfTransitions())
{

	typedef typename Automaton<Y>::Weight Weight2;


	typename Automaton<Weight2>::TransitionMatrix::const_iterator amit = a.trans_matrix.begin();

	for(; amit!= a.trans_matrix.end(); ++amit) {
		int letter = amit->first;
		typename Automaton<Weight2>::TransitionNZ::const_iterator ait = a.trans_nz.find(letter);

		if(ait==a.trans_nz.end()) continue;

		trans_nz[letter] = ait->second;

		trans_matrix[letter] = amit->second;
	}

	accepting_actions.insert(a.accepting_actions.begin(),a.accepting_actions.end());
}


template<typename X>
Automaton<X>::Automaton(const Automaton& a, const Automaton& b) :
	initial(a.nr_of_states + b.nr_of_states,1),
	accepting(a.nr_of_states + b.nr_of_states,1),
	nr_of_states (a.nr_of_states + b.nr_of_states),
	nr_of_transitions (a.nr_of_transitions + b.nr_of_transitions)
{

	VectorRange initial_a (initial,boost::numeric::ublas::range(0,a.nr_of_states));
	VectorRange initial_b (initial,boost::numeric::ublas::range(a.nr_of_states,nr_of_states));
	initial_a = a.initial;
	initial_b = b.initial;

	VectorRange accepting_a (accepting,boost::numeric::ublas::range(0,a.nr_of_states));
	VectorRange accepting_b (accepting,boost::numeric::ublas::range(a.nr_of_states,nr_of_states));
	accepting_a = a.accepting;
	accepting_b = - b.accepting;

	std::unordered_set<int> letters;


	for(typename TransitionMatrix::const_iterator amit = a.trans_matrix.begin() ; amit != a.trans_matrix.end(); ++ amit) {
	  letters.insert(amit->first);
	}

	for(typename TransitionMatrix::const_iterator bmit = b.trans_matrix.begin() ; bmit != b.trans_matrix.end(); ++ bmit) {
	  letters.insert(bmit->first);
	}

	for(std::unordered_set<int>::const_iterator it = letters.begin(); it != letters.end(); ++ it) {

	  int i = *it;

	  std::cout << ">>>>>>> Letter " << alphabet[i] << std::endl;

		std::vector<unsigned>& nz(trans_nz[i]);

		unsigned a_transnz_size = 0;
		TransitionNZ::const_iterator ait = a.trans_nz.find(i);

		if(ait!=a.trans_nz.end()) {
			nz.insert(nz.end(), ait->second.begin(), ait->second.end());
			a_transnz_size = ait->second.size();
		}

		TransitionNZ::const_iterator bit = b.trans_nz.find(i);
		unsigned b_transnz_size = 0;

		bit = b.trans_nz.find(i);
		if(bit!=b.trans_nz.end()) {
			nz.insert(nz.end(), bit->second.begin(), bit->second.end());
			b_transnz_size = bit->second.size();
		}

		if( ait!=a.trans_nz.end() && bit!=b.trans_nz.end() ) {
			for(unsigned j=ait->second.size(); j<nz.size(); ++j)
				nz[j] += a.nr_of_states;
		}


		trans_matrix[i].resize(trans_nz.size(),nr_of_states,false);

		Matrix& M(trans_matrix[i]);

		M.resize( a_transnz_size + b_transnz_size, nr_of_states, false );

		typename TransitionMatrix::const_iterator amit = a.trans_matrix.find(i);
		if(amit != a.trans_matrix.end()) {
			MatrixRange M_a (M,boost::numeric::ublas::range(0,amit->second.size1()),boost::numeric::ublas::range(0,amit->second.size2()));
			M_a = amit->second;
		}

		typename TransitionMatrix::const_iterator bmit = b.trans_matrix.find(i);
		if(bmit != b.trans_matrix.end()) {
			MatrixRange M_b (M,boost::numeric::ublas::range(amit->second.size1(), M.size1()), boost::numeric::ublas::range(amit->second.size2(),M.size2()));
			M_b = bmit->second;
		}
	}

	enabled.reserve(a.enabled.size() + b.enabled.size());
	enabled.insert(enabled.end(),a.enabled.begin(), a.enabled.end());
	enabled.insert(enabled.end(),b.enabled.begin(), b.enabled.end());

	incoming.reserve(a.incoming.size() + b.incoming.size());
	incoming.insert(incoming.end(),a.incoming.begin(), a.incoming.end());
	incoming.insert(incoming.end(),b.incoming.begin(), b.incoming.end());

	accepting_actions.insert(a.accepting_actions.begin(),a.accepting_actions.end());
	accepting_actions.insert(b.accepting_actions.begin(),b.accepting_actions.end());
}


// count the non-zero rows

template<typename X>
inline
  long Automaton<X>::getMem() const {
	  long result = 0;

	  // memory for the matrics

	  for(typename std::unordered_map<int,Matrix>::const_iterator mit = trans_matrix.begin(); mit!= trans_matrix.end(); ++mit) {
	    const Matrix& M = mit -> second;
	    result += mem(M);
	  }
	  result += sizeof(trans_matrix);

	  // memory for the vectors
	  result += mem(accepting);
	  result += mem(initial);


	  for(unsigned i = 0 ; i<enabled.size(); ++i) {
	    result += sizeof(int) * enabled[i].size();
	    result += sizeof(enabled[i]);
	  }

	  for(unsigned i = 0 ; i<incoming.size(); ++i) {
	    result += sizeof(int) * incoming[i].size();
	    result += sizeof(incoming[i]);
	  }

	  result += sizeof(Automaton);

	  return result;
  }



template<typename X>
inline
void Automaton<X>::computeMatrices() {

  TransitionComp tc;

  // sort the transitions according to their actions
  std::sort(trans.begin(),trans.end(),tc);

  for(typename std::vector<Transition>::const_iterator it=trans.begin(); it!=trans.end(); ++it) {
    const Transition& t (*it);

	// check if row is already accounted for
	if(trans_nz[t.sigma].size() > 0) {
		unsigned previous = trans_nz[t.sigma].size()-1;
		// if so do not add it
		if(trans_nz[t.sigma][previous] == t.source)
			continue;
    	}

	trans_nz[t.sigma].push_back(t.source);
  }

  for(unsigned i=0; i<alphabet.size();++i) {

	if(trans_nz[i].size() > 0)
	trans_matrix[i].resize(trans_nz[i].size(),nr_of_states,false);
  }

  int prev_sigma = -1;
  unsigned prev_source = 0;
  unsigned row_index = 0;

  for(typename std::vector<Transition>::const_iterator it=trans.begin();
	  it!=trans.end(); ++it)  {
	const Transition& t (*it);
	int target = t.target;
	X prob = t.prob;

	if(prev_sigma != t.sigma) {
		row_index = 0;
		prev_source = t.source;
	}

	if( prev_source != t.source )
		++row_index;

	enabled[t.source].insert(t.sigma);
	incoming[target].insert(t.sigma);

    	trans_matrix[t.sigma](row_index,target) += prob; // sum up multiple transitions with the same label and target


	prev_source = t.source;
	prev_sigma = t.sigma;
  }

  trans.clear(); 
}

template<typename X>
inline
unsigned Automaton<X>::getNrOfNonzeros() {

  unsigned result = 0;
  for(typename std::unordered_map<int,Matrix>::const_iterator mit = trans_matrix.begin(); mit!= trans_matrix.end(); ++mit) {
    const Matrix& M = mit -> second;
    for(typename Matrix::const_iterator1 it1 = M.begin1(); it1!=M.end1(); ++it1) {
      for(typename Matrix::const_iterator2 it2 = it1.begin(); it2!=it1.end(); ++it2) {
		if(*it2 != 0) ++result;
      }
    }
  }
  return result;
}

template<typename X>
inline
void Automaton<X>::addTransition(unsigned int source, unsigned int target, std::string sigma, X prob){

  int index;


  std::unordered_map<std::string,int>::const_iterator it(alphabet_table.find(sigma));
  if(it!=alphabet_table.end()) {
    index = it -> second;
  } else {
    alphabet_table[ sigma ] = index = alphabet.size();
    alphabet.push_back(sigma);
  }

  addTransition(source, target, index, prob);
}

template<typename X>
inline
void Automaton<X>::addTransition(unsigned int source, unsigned int target, int sigma, X prob){

  int index = sigma;

  static Transition t;

  assert(source < nr_of_states);
  assert(target < nr_of_states);

  t.source = source;
  t.target = target;
  t.sigma = index;

  if(prob==0)
    return;

  if(std::numeric_limits<X>::infinity() == prob ) {
    std::cerr << "Automaton::addTransition source state " << source << " target " << target << " invalid transition probability " << prob << std::endl;
    exit(1);
  }

  t.prob = prob;


  if(target == accepting_state)
    accepting_actions.insert(index);

  trans.push_back(t);


  ++nr_of_transitions;

}

template<typename X>
inline
bool Automaton<X>::checkRangeOfStateIndex(int index) const {
	bool result = true;
	if(index < 0) {
		std::cout << "State index negative : " << index << std::endl;
		result = false;
	} else if( (unsigned) index >= nr_of_states) {
		std::cout << "State index out of range : " << index << " > " << nr_of_states << std::endl;
		result = false;
	}

	return result;
}



/**
  * reads automaton from file
  * the format is as follows:
  *
  * <initial> <nr_of_states>
  * <accepting state 1> [(<val>,<num>,<deno>)]*
  * ...
  * <accepting state k> [(<val>,<num>,<deno>)]*
  * -1
  * <source> <target> <letter> <num> <deno>
  * ...
  * -1
  */
template<typename X>
inline
int Automaton<X>::parseTransition(char* input) {
  int source, target;
  long num, deno;
  double prob1;
  X prob;
  static char sigma[1000];
  if(sscanf(input,"%d %d %s %lu %lu",&source,&target,sigma,&num,&deno) == 5) {
	if(deno == 0)
	  	std::cout << "parseTransition: division by zero" << std::endl;

    prob = ((X)num) / ((X)deno);
  }
  else if(sscanf(input,"%d %d %s %lg",&source,&target,sigma,&prob1) == 4) {
	prob = ((X)prob1);
  } else {
    std::string s(input);
    throw "unable to read input \'" + s + "\'";
  }

  try {
	  checkRangeOfStateIndex(source);
	  checkRangeOfStateIndex(target);

	  addTransition(source,target,sigma,prob);
  } catch (std::string msg) {
  	  std::cerr << "Automaton::readTransition: failed to read " << input << std::endl;
	  throw msg;
  }
  return 0;
}

template<typename X>
inline
int Automaton<X>::parseDistribution(char* input, int accepting_state) {

  int state, value;
  long num, deno;
  double prob1;
  X prob;

  // ignore white space
  while(*input == ' ') { ++input; }
  
  if(sscanf(input,"%d",&state)!=1)
	  return -1;

  while(*input != '\0') {

    /* skip to next whitespace */
    while( *input !=' ' && *input != '\0') {
      ++input;
    }

    /* skip next white space */
    while( *input ==' ' && *input != '\0') {
      ++input;
    }



    if(sscanf(input,"(%d,%lu,%lu)",&value,&num,&deno) == 3) {
      if(deno == 0)
    	  std::cout << "parseDistribution: division by zero" << std::endl;
      prob = ((X)num) / ((X)deno);
    }
    else if(sscanf(input,"(%d,%lg)",&value,&prob1) == 2) {
      prob = ((X)prob1);
    } else {
      return -1;
    }

    if(value == -1)
      break;
    else {
      static char buf[1000];
      sprintf(buf,"%d",value);
      addTransition(state,accepting_state,buf,prob);
    }
  }
  return 0;
}

template<typename X>
inline
int Automaton<X>::readVector(char* input, Automaton<X>::Vector& vec) const {

  int value;
  long num, deno;
  double prob;

    /* skip to next whitespace */
    while( *input ==' ' && *input != '\0') {
      ++input;
    }

  
  while(*input != '\0') {

    if(sscanf(input,"%lu/%lu:%d",&num,&deno,&value) == 3) {
       vec(value)= X( num ) / X( deno );
    }
    // expecting <real number> : <state number>
    else if(sscanf(input,"%lg",&prob) == 1) {
    	    
     /* skip non-white space */
     while( (*input!= ':') && (*input != '\0')) {
       ++input;
      }	    
    	    
      if(*input==':')
        ++input;

      /* skip whitespace */
      while( *input ==' ' && *input != '\0') {
        ++input;
      }

       if(sscanf(input,"%d",&value)==1)
       {
         vec(value) = X(prob);
       } else
       {
       	       std::string s(input);
       	       throw "Vector " + s + " could not match ... : <state>";
       }
    } else {
       std::string s(input);
       throw "Vector " + s + " could not be parsed";
    }
    /* skip non-white-space characters*/
    while( !isspace(*input) && (*input != '\0')) {
       ++input;
    }

    /* skip to next whitespace */
    while( *input ==' ' && *input != '\0') {
      ++input;
    }
  }
  return 0;
}

template<typename X>
inline
int Automaton<X>::writeVector(const Vector& v, std::ofstream& f) {
  for(typename Vector::const_iterator vit = v.begin(); vit!= v.end(); ++vit) {
	if(*vit != 0)
    f << *vit << ":"<<vit.index()<< " ";
  }
  return 0;
}


template<typename X>
inline
int Automaton<X>::writeTransitions(std::ofstream& f) {
 for(typename std::unordered_map<int,Matrix>::const_iterator mit = trans_matrix.begin(); mit!= trans_matrix.end(); ++mit) {
    int action = mit -> first;
    const Matrix& M = mit -> second;
    const std::string& letter( alphabet[action] );

    const std::unordered_map<int,std::vector<unsigned> >::const_iterator nzit (trans_nz.find(action));

    if(nzit == trans_nz.end() ) return -1;
    const std::vector<unsigned>& nz = nzit -> second;

    // iterate over the transitions
    for(typename Matrix::const_iterator1 it1 = M.begin1(); it1!= M.end1(); ++it1) {
      for(typename Matrix::const_iterator2 it2 = it1.begin(); it2!=it1.end(); ++it2) {
    	  f << nz[it2.index1()] << " " << it2.index2() << " " << letter << " " << *it2 << "\n";
      }
    }
  }
 return 0;
}

template<typename X>
inline
int Automaton<X>::writeMatlabVector(const Vector& v, std::ofstream& f) {
  f << "[";
  for(typename Vector::const_iterator vit = v.begin(); vit!= v.end(); ++vit) {
    if(vit != v.begin())
	f << " ";
    f << *vit;
  }
  f << "]";
  return 0;
}

template<typename X>
inline
int Automaton<X>::writeMatlabMatrix(std::ofstream& f) {



 for(typename std::unordered_map<int,Matrix>::const_iterator mit = trans_matrix.begin(); mit!= trans_matrix.end(); ++mit) {
    int action = mit -> first;
    const Matrix& M = mit -> second;
    const std::string& letter( alphabet[action] );

    const std::unordered_map<int,std::vector<unsigned> >::const_iterator nzit (trans_nz.find(action));

    if(nzit == trans_nz.end() ) return -1;

    f << "A" << letter << "=[";


    const std::vector<unsigned>& nz = nzit->second;

    int nzctr = 0;
    typename Matrix::const_iterator1 it1 = M.begin1();    

    for(unsigned i = 0; i<nr_of_states; ++i) {

      if(i>0)
	f << ";\n";

      if(nz[nzctr] == i) { // found a non-zero row

	for(unsigned j=0; j<M.size2(); ++j) {
    	  f << M(it1.index1(),j);
    	  if(j + 1 < M.size2()) {
    		  f << " ";
	  }
	}

	++nzctr;
	++it1;
      } else { // this is a zero row 
	// print a bunch of zeros
	for(unsigned j = 0; j<nr_of_states; ++j) {
	  f << "0" ;
	  if( j + 1 < nr_of_states)
	    f << " ";
	}
      }
   }
   f << "];";
  }
 return 0;
}


template<typename X>
inline
int Automaton<X>::write(const char* filename) {
  std::ofstream f (filename);
  unsigned n (getNrOfStates());

  // set the precision for output of floating point numbers
  f.setf(std::ios_base::floatfield);
  f.precision(std::numeric_limits<double>::digits10 );

  /***** number of states *****/
  f << n << "\n";

  f << "initial\n";
  writeVector(initial,f);
  f << "\n";

  f << "accepting\n";
  writeVector(accepting,f);
  f << "\n";

  f << "trans\n";
  writeTransitions(f);

  return 0;
}

template<typename X>
inline
int Automaton<X>::writeApex(const char* filename) {
  std::ofstream f (filename);
  unsigned initial_state (0);
  unsigned n (getNrOfStates());

  // set the precision for output of floating point numbers
  f.setf(std::ios_base::floatfield);
  f.precision(std::numeric_limits<double>::digits10);

  f << initial_state << " " << n << "\n";

  /***** accepting states and their distributions *****/


  f << "-1" << std::endl; // end the accepting-state section

  /***** transitions *****/
  writeTransitions(f);

  f << "-1" << std::endl; // end the transition section

  return 0;
}

template<typename X>
inline
int Automaton<X>::writeMatlab(const char* filename) {
  std::ofstream f (filename);

  // set the precision for output of floating point numbers
  f.setf(std::ios_base::floatfield);
  f.precision(std::numeric_limits<double>::digits10);

  writeMatlabVector(initial,f); f << "\n";
  writeMatlabVector(accepting,f); f << "\n";


  /***** transitions *****/
  writeMatlabMatrix(f);

  return 0;
}

template<typename X>
inline
int Automaton<X>::readApex(const char* filename) {
   int result = 0;

   unsigned line = 0;

   try {
	  std::ifstream f (filename);

	  if(!f) {

		printf("Wrong command-line option: file \"%s\" does not exist\n",filename);
		return -1;
	  }

	  char buffer[buffer_length];

	  f.getline(buffer,buffer_length);
	  ++line;

	  int i,n;
	  if(sscanf(buffer,"%d %d",&i,&n) != 2) return -1;

	  if( i >= n ) {
		  throw ("Line " + toString(line) + ": invalid initial state " + toString(i));
	  }

	  setNrOfStates(n + 1); // !!! why +1 : because this is the accepting state
	  setInitial(i);
	  setAccepting(n);

	  while(true) {

		f.getline(buffer,buffer_length);
		++line;
		if(strcmp(buffer,"-1") == 0)
		  break;
		result = parseDistribution(buffer,n);
		if(result)
			throw "Line " + toString(line) + " : could not parse distribution \'" + buffer + "\'";
	  }
	  do {
		f.getline(buffer,buffer_length);
		++line;
		if(strcmp(buffer,"-1") == 0)
		  break;
		try {
			result = parseTransition(buffer);
		} catch (std::string msg) {
			throw "Line " + toString(line) + " : " + msg;
		}
	  } while (strcmp(buffer,"-1"));

          if(global_setting::verbosity >= 1)
	  std::cout<<"   # states: " << n << " # trans: " << trans.size() << " # letter " << alphabet.size()<< std::endl;

	   f.close();

   } catch(std::string msg) {
	   std::cout <<msg<<std::endl;
	   result = -1;
   }

  return result;
}

template<typename X>
inline
int Automaton<X>::read(const char* filename) {
	int result = 0;

  unsigned line = 1;

  try {
	  std::ifstream f (filename);

	  if(!f) {

		printf("Wrong command-line option: file \"%s\" does not exist\n",filename);
		return -1;
	  }

	  char buffer[buffer_length];

	  f.getline(buffer,buffer_length);

	  ++line;

	  int n;
	  if(sscanf(buffer,"%d",&n) != 1) throw "Line " + toString(line) + " : Number of states not specified";

	  setNrOfStates(n);

	  f.getline(buffer,buffer_length);

	  ++line;

	  if(strcmp(buffer,"initial") == 0) {
		  f.getline(buffer,buffer_length);


		  try {
			result = readVector(buffer,initial);
			++line;
		  } catch (std::string msg) {
			  throw "Line " + toString(line) + " : initial " + msg;
		  }
	  } else throw "Line " + toString(line) + " : Keyword \'initial\' is missing";;

	  f.getline(buffer,buffer_length);
	  ++line;
	  if(strcmp(buffer,"accepting") == 0) {
		  f.getline(buffer,buffer_length);

		  try {
			result = readVector(buffer,accepting);
			++line;
		  } catch (std::string msg) {
			  throw "Line " + toString(line) + " : accepting " + msg;
		  }

	  } else throw "Line " + toString(line) + " : Keyword \'accepting\' is missing";

	  if(result)
		  return result;

	  f.getline(buffer,buffer_length);
	  ++line;

	  if(strcmp(buffer,"trans") != 0) return -1;
	  while(!f.eof()) {
		f.getline(buffer,buffer_length);
		if(strcmp(buffer,"") == 0) {
			break;
		}

		try {
			result = parseTransition(buffer);
		} catch (std::string msg) {
			throw "Line " + toString(line) + " : " + msg;
		}

		++line;
	  }
	  f.close();
	  if(global_setting::verbosity >= 1)
	  	std::cout<<"   # states: " << n << " # trans: " << trans.size() << " # letter " << alphabet.size()<< std::endl;
  } catch (std::string msg) {
	  if(global_setting::verbosity >= 1)
		  std::cout << msg << std::endl;
	  result = -1;
  } catch (char* msg) {
	  if(global_setting::verbosity >= 1)
	  	std::cout << msg << std::endl;
	  result = -1;
  }

  return result;
}

template<typename X>
inline
X Automaton<X>::operator()(const Word& w) const {

  Vector v(initial);

  for(Word::const_iterator it=w.begin(); it!=w.end();++it) {

    // get the corresponding action
    int sigma(*it);

    // fetch the matrix
    typename std::unordered_map<int,Matrix>::const_iterator mit = trans_matrix.find(sigma);

    if( mit == trans_matrix.end()) {
      std::cout << "Automaton::operator(): trace contains non-existant action " << alphabet[sigma] << std::endl;
      return 0.0;
    }

  	Vector vnew (v.size());
	vnew.clear();

	// fetch the nz entries
	std::unordered_map<int,std::vector<unsigned> >::const_iterator nit = trans_nz.find(sigma);

    if( nit == trans_nz.end()) {
      std::cout << "Automaton::operator(): unexpected error, missing non-zero entry " << std::endl;

      return 0.0;
    }

    row_prod(vnew,v,mit->second,nit->second);
    v = vnew;
  }

  // multiply with accepting state vector
  return prec_inner_prod(v,accepting);
}


inline
void printAiSeeNode(std::ostream &stream, int id, const std::string& label,
		    const std::string& infobox, const std::string& infobox2, int width,
		    int height,
		    const std::string& shape,
		    const std::string& fill,
		    const std::string textcolor = "",
		    int borderwidth = -1) {
  stream << "node: { " << "\n" << "  title: \"" << id << "\"\n"
	 << "  width: " << width << "\n" << "  height: " << height << "\n"
	 << "  shape: " << shape << "\n";
  if ("" != infobox) {
    stream << "  info1: \"" << infobox << "\"\n";
  }
  if ("" != infobox2) {
    stream << "  info2: \"" << infobox2 << "\"\n";
  }
  if ("" != fill) {
    stream << "  color: " << fill << "\n";
  }
  if ("" != textcolor) {
    stream << "  textcolor: " << textcolor << "\n";
  }


  if(-1 != borderwidth) {
    stream << "  borderwidth: " << borderwidth << "\n";
  }

  stream << "  label: \"" << label << "\"\n" << "}\n";
}

/**
 * Prints graph edge with given attributes.
 */
inline
void printAiSeeEdge(std::ostream &stream, int source, int target,
		    const std::string& label, int width,
		    // "line", "dashed", "dotted"
		    const std::string& style,
		    //"delta", "standard", "diamond", "short", "white_delta", "white_diamond", or "none"
		    const std::string& sourceArrow,
		    //"delta", "standard", "diamond", "short", "white_delta", "white_diamond", or "none"
		    const std::string& targetArrow, const std::string& fill) {
  stream << "edge: {\n " << "  sourcename: \"" << source << "\"\n"
	 << "  targetname: \"" << target << "\"\n" << "  thickness: "
	 << width << "\n";
  if ("" != style) {
    stream << "  style: " << style << "\n";
  }
  if ("" != fill) {
    stream << "  color: " << fill << "\n";
  }
  stream << "  backarrowstyle: " << sourceArrow << "\n" << "  arrowstyle: "
	 << targetArrow << "\n" << "  label: \"" << label << "\"\n"
	 << "  class: 1}\n";
}

template<typename X>
inline
void Automaton<X>::aiSee(std::ostream& stream) const {

  stream << "graph: {" << std::endl << std::endl
	 << "colorentry 1: 153 255 0" << std::endl
	 << "colorentry 2: 255 255 0" << std::endl
	 << "colorentry 3: 0 0 0" << std::endl // black
	 << "colorentry 4: 255 255 221" << std::endl
	 << "colorentry 5: 128 0 128" << std::endl // purple
	 << "colorentry 6: 0 0 255" << std::endl // blue
	 << "colorentry 7: 200 0 0" << std::endl << std::endl // red
	 << "colorentry 8: 200 200 200" << std::endl // lightgray
	 << "colorentry 9: 255 255 255" << std::endl; // white

  // shapes
  const std::string box("box");
  const std::string rhomb("rhomb");
  const std::string hexagon("hexagon");
  const std::string triangle("triangle");
  const std::string ellipse("ellipse");

  // colors
  const std::string green("1");
  const std::string yellow("2");
  const std::string black("3");
  const std::string gray("4");
  const std::string purple("5");
  const std::string blue("6");
  const std::string red("7");
  const std::string lightgray("8");
  const std::string white("9");

  int width, height;

  std::string label;
  std::string sourceArrow = "none";
  std::string targetArrow = "solid";
  std::string style;
  std::string shape;
  std::string fill = black;

  std::string infobox1, infobox2;

  // iterate through the states
  for (unsigned i = 0; i < nr_of_states; ++i) {
    std::string fillcolor;
    std::string textcolor(black);
    int borderwidth = -1;

    label = toString(i);

    shape = ellipse;
    if(accepting(i) > 0) {
      fillcolor = red;
      textcolor = white;
      width = -1;
      height = -1;
    } else if(initial(i) > 0) {
      fillcolor = green;
      width = -1;
      height = -1;
    } else {
      fillcolor = gray;
      width = -1;
      height = -1;
    }

    printAiSeeNode(stream, i, label,
		   infobox1, // infobox1
		   infobox2, // infobox2
		   width, // width
		   height, // height
		   shape,
		   fillcolor,
		   textcolor,
		   borderwidth
		   );
  }


  // iterate through the actions
  for(typename std::unordered_map<int,Matrix>::const_iterator mit = trans_matrix.begin(); mit!= trans_matrix.end(); ++mit) {
    int action = mit -> first;
    const Matrix& M = mit -> second;
    const std::string& letter( alphabet[action] );


    const std::unordered_map<int,std::vector<unsigned> >::const_iterator nzit (trans_nz.find(action));

    if(nzit == trans_nz.end() ) return;

    const std::vector<unsigned>& nz = nzit -> second;

    // iterate over the transitions
    for(typename Matrix::const_iterator1 it1 = M.begin1(); it1!= M.end1(); ++it1) {
      for(typename Matrix::const_iterator2 it2 = it1.begin(); it2!=it1.end(); ++it2) {
	label = toString(*it2) + "," + letter;
	width = 1;
	fill = "3";

	printAiSeeEdge(stream, nz[it2.index1()], it2.index2(), label, width, style,
		       sourceArrow, targetArrow, fill);
      }
    }
  }



  stream << "}" << std::endl;
}


#endif
