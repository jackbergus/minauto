#ifndef MUTATE
#define MUTATE

template<class Automaton>
void mutateAutomaton(Automaton& a, unsigned steps = 1);


// IMPLEMENTATION


template<class Automaton>
void mutateAutomaton(Automaton& a, unsigned steps) {

  typedef typename Automaton::Weight Scalar;
  typedef typename Automaton::TransitionMatrix::iterator TransitionMatrixIterator;
  typedef typename Automaton::Matrix Matrix;

  typedef typename Matrix::iterator1 MatrixIterator1;
  typedef typename Matrix::iterator2 MatrixIterator2;

  // select an action at random
  unsigned nr_of_actions (a.trans_matrix.size()); // how many actions are there


  for(unsigned step = 0 ; step < steps ; ++step) {

    unsigned rand_action (rand() % nr_of_actions); // pick action

    // pick the corresponding matrix
    TransitionMatrixIterator it=a.trans_matrix.begin();
    for(unsigned counter=0;it != a.trans_matrix.end() && counter < rand_action; ++it,++counter);
    Matrix& M(it->second);

    // select a non-zero of the corresponding matrix at random
    unsigned nnz = 0;
    for(MatrixIterator1 it1 = M.begin1(); it1!= M.end1(); ++it1) {
      for(MatrixIterator2 it2 = it1.begin(); it2!= it1.end(); ++it2) {
    	  Scalar prob(*it2);
    	  if(prob == 0.0) continue;
			  ++nnz;
		  }
		}

    unsigned rand_nz_choice (rand());
    std::cout << "rand choice " << rand_nz_choice << " choice " << rand_nz_choice % nnz << " nnz " << nnz << " " << std::endl;
    unsigned rand_nz ( rand_nz_choice % nnz );

    // modify the non-zero entry "rand_nz" to a random probability "rand_prob"

    double rand_prob (getRandomDouble());

    unsigned i = 0,j =0;

    bool done = false;
    unsigned counter = 0;
    for(MatrixIterator1 it1 = M.begin1(); !done && it1!= M.end1(); ++it1) {
      i = it1.index1();

      for(MatrixIterator2 it2 = it1.begin(); !done && it2!= it1.end(); ++it2) {
	j = it2.index2();

	Scalar& prob(*it2);
	if(prob == 0.0) continue;
	if(counter == rand_nz) {
	  prob = rand_prob;
	  done = true;
	  break;
	}
	++counter;
      }
    }
    std::cout << "Mutated matrix of action "<<alphabet[it->first] << "("<<i<<","<<j<<") = "<<rand_prob <<std::endl;
  }


}


#endif
