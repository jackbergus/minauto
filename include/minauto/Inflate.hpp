#ifndef INFLATE
#define INFLATE


template<class Automaton>
void forwardInflate(Automaton& a, unsigned steps);


template<class Automaton>
void backwardInflate(Automaton& a, unsigned steps);


// IMPLEMENTATION


template<class Vector>
void selectRandomVector(Vector& v, Vector& acc) {
  // select random entries
  for(unsigned i=0; i<v.size(); ++i) {
    v[i] = acc(i) != 0.0 ? 0.0 : getRandomDouble(10);
  }

  // enfore that the vector becomes a distribution, i.e., "sums up to one"
  v /= norm_2(v);

}

/*
  scale the vectors col such that M - outer_prod(col,row)
  is a (sub)stochastic matrix
 */
template<class Matrix, class Vector>
void ensurePositiveWeights(Matrix& M, Vector& delta, Vector& ssigma) {

	/*
  unsigned dim = M.size1();



  for(Matrix::const_iterator1 it=M.begin1(); it!=M.end1() ; ++it) {
    double ssigma_row ( ssigma(it.index1()));

	unsigned counter = 0;

    for(Matrix::const_iterator2 it2 = it.begin(); it2!=it.end(); ++it2) {
      double delta_col (delta (it2.index2()));
      double matrix_entry ( *it2 );
	  ++counter;
      if(matrix_entry == 0.0 && row_entry > 0.0 ) {
		col(it2.index2()) = 0.0;
      } else {
		assert(false);
		double p(row_entry * col_entry);
		if ( p > matrix_entry) {
		  col(it2.index2()) *= p / matrix_entry;
		}
      }
    }

	zero |= counter < M.size2();
  }


  for(SparseVector::iterator it = col.begin(); it!=col.end(); ++it) {
    *it = 0;
  }
  */
}

template<class Matrix, class Vector>
void ensurePositiveWeights2(Matrix& M, Vector& col, Vector& row) {

  /*
  for(SparseVector::iterator it = col.begin(); it!=col.end(); ++it) {
    *it = 0;
  }

  unsigned dim = M.size1();
  for(unsigned i=0; i<dim ; ++i) {
    double row_entry ( row(i));

    for(unsigned j=0; j<dim ; ++j) {
      double col_entry (col (j));
      double matrix_entry ( M(i,j) );
      if(matrix_entry == 0.0 && row_entry > 0.0 ) {
	row(j) = 0.0;
      } else {
	double p(row_entry * col_entry);
	if ( p > matrix_entry) {
	  row(j) *= p / matrix_entry;
	}
      }
    }
  }
*/
}


template<class Automaton>
inline
void forwardInflate(Automaton& a, unsigned steps) {

  typedef typename Automaton::Vector Vector;
  typedef typename Automaton::Matrix Matrix;
  typedef typename Automaton::TransitionMatrix::iterator TransitionMatrixIterator;

  unsigned n(a.getNrOfStates());

  std::cout << "Forward inflate: nr of non-zeros (before) " << a.getNrOfNonzeros() << std::endl;

  // unfold into uncompressed form with potential zero rows
  for(TransitionMatrixIterator it=a.trans_matrix.begin();
	it != a.trans_matrix.end(); ++it) {
    int sigma (it->first);
	Matrix M_sigma ( it -> second );

	unfold(M_sigma, a.trans_nz[sigma], it->second);
  }


  for(unsigned i=0; i<steps; ++i) {

    Vector delta (n);
    selectRandomVector<Vector>(delta,a.accepting);



    /*
                                                          column n+1
                       -----------------------------------------------
                       |   M(sigma) - s(sigma) x delta  |   s(sigma) |
      M'(sigma) =      |________________________________|            |
                       |   delta x M(sigma)                    0     |             row (n+1)
		       -----------------------------------------------

    */
    for(TransitionMatrixIterator it=a.trans_matrix.begin();
	it != a.trans_matrix.end(); ++it) {
      int sigma (it->first);


      Matrix& M_sigma ( it -> second );
      Vector  s_sigma (n) ;

      Vector delta_times_M_sigma(n);

      if(a.accepting_actions.find(sigma) == a.accepting_actions.end()) {
		selectRandomVector(s_sigma, a.accepting);

		ensurePositiveWeights(M_sigma, delta, s_sigma);

		/* compute delta x M(sigma) */
		delta_times_M_sigma = prec_prod(delta,M_sigma);

		/* compute M(sigma) - s(sigma) x delta */
		M_sigma.minus_assign(outer_prod(s_sigma,delta));
      }

      assert(M_sigma.size1() == n && M_sigma.size2() == n);


      /* increase the dimension of M_sigma */

      Matrix tmp(M_sigma);

      M_sigma.resize(n+1,n+1,false);

      project(M_sigma,boost::numeric::ublas::range(0,n),boost::numeric::ublas::range(0,n)).assign(tmp);

	  assert(M_sigma.size1() == n+1);


      //std::cout << "M(sigma) (nnz) "<<  M_sigma.nnz() << " " << M_sigma.size1() << " x " << M_sigma.size2()<< " delta x M(sigma) nnz " << delta_times_M_sigma.nnz() << " delta (nnz) "<< delta.nnz()<< " s(sigma) " << ssigma.nnz() << std::endl;


      /* add last row */

      delta_times_M_sigma.resize(n + 1,true); // make sure that the row has the right dimension
      delta_times_M_sigma(n) = 0;
      row(M_sigma,  n) = delta_times_M_sigma;

      /* add last column */
      s_sigma.resize(n + 1,true);
      s_sigma(n) = 0;
      column(M_sigma,  n) = s_sigma;

      //std::cout << "M(sigma) (nnz) "<<  M_sigma.nnz() <<std::endl;

      /* walk through last row to update incoming */

      const Vector& last_row = row(M_sigma,n);
      for(typename Vector::const_iterator it=last_row.begin(); it!=last_row.end(); ++it) {
		if(*it == 0.0) continue;
		a.incoming[it.index()].insert(sigma);
      }

      /* walk through last column to update outgoing */
      const Vector& last_column = column(M_sigma,n);

      for(typename Vector::const_iterator it=last_column.begin(); it!=last_column.end(); ++it) {
		if(*it == 0.0) continue;
		a.enabled[it.index()].insert(sigma);
      	}
      }

      /* increase dimension of initial and accepting state vector */
      delta.resize(n+1,true);
	  delta(n) = 0.0;

      a.setNrOfStates(n+1);

      a.accepting(n) = prec_inner_prod(delta,a.accepting);
      n = n + 1;
  }


  // fold into compressed form (without zero rows)
  for(TransitionMatrixIterator it=a.trans_matrix.begin();
	it != a.trans_matrix.end(); ++it) {
    int sigma (it->first);
	Matrix M_sigma ( it -> second );

	fold(M_sigma, it->second, a.trans_nz[sigma]);
  }


  std::cout << "Forward Inflate: nr of non-zeros (after)" << a.getNrOfNonzeros() << std::endl;
}

template<class Automaton>
inline
void backwardInflate(Automaton& a, unsigned steps) {
  typedef typename Automaton::Vector Vector;
  typedef typename Automaton::Matrix Matrix;
  typedef typename Automaton::TransitionMatrix::iterator TransitionMatrixIterator;

  unsigned n(a.getNrOfStates());


  // unfold into uncompressed form with potential zero rows
  for(TransitionMatrixIterator it=a.trans_matrix.begin();
	it != a.trans_matrix.end(); ++it) {
    int sigma (it->first);
	Matrix M_sigma ( it -> second );

	unfold(M_sigma, a.trans_nz[sigma], it->second);
  }

  std::cout << "Backward inflate: nr of non-zeros (before) " << a.getNrOfNonzeros() << std::endl;


  for(unsigned i=0; i<steps; ++i) {

    Vector delta (n);
    selectRandomVector<Vector>(delta,a.accepting);


    //std::cout << std::setprecision(10);
    //std::cout << "a.accepting " << a.accepting << std::endl;
    //std::cout << "delta " << delta << std::endl;


    /*
                                                          column n+1
                       ----------------------------------------------------
                       |   M(sigma) - delta x s(sigma)  |   M(sigma)delta |
      M'(sigma) =      |________________________________|                 |
                       |   s(sigma)                                 0     |             row (n+1)
		       ----------------------------------------------------

    */
    for(TransitionMatrixIterator it=a.trans_matrix.begin();
	it != a.trans_matrix.end(); ++it) {
      int sigma (it->first);

      Matrix& M_sigma ( it -> second );
      Vector  s_sigma (n) ;

      Vector M_sigma_times_delta(n);

      if(a.accepting_actions.find(sigma) == a.accepting_actions.end()) {
		selectRandomVector<Vector>(s_sigma, a.accepting);
		ensurePositiveWeights2(M_sigma, delta, s_sigma);

		/* compute M(sigma) x delta */
		M_sigma_times_delta = prod(M_sigma,delta);

		/* compute M(sigma) - delta x sigma */
		M_sigma.minus_assign(outer_prod(delta,s_sigma));

		//std::cout << "s_sigma " << s_sigma << std::endl;
		//std::cout << "outer_prod(delta,s_sigma) " << outer_prod(delta,s_sigma) << std::endl;
		//std::cout << "M_sigma" << M_sigma << std::endl;

      }

      assert(M_sigma.size1() == n && M_sigma.size2() == n);


      /* increase the dimension of M_sigma */

      Matrix tmp(M_sigma);

      M_sigma.resize(n+1,n+1,false);

      //std::cout << " " << n + 1 << " " <<
      project(M_sigma,boost::numeric::ublas::range(0,n),boost::numeric::ublas::range(0,n)).assign(tmp);

      Vector ssigma(s_sigma);
      //std::cout << "M(sigma) (nnz) "<<  M_sigma.nnz() << " " << M_sigma.size1() << " x " << M_sigma.size2()<< " delta x M(sigma) nnz " << delta_times_M_sigma.nnz() << " delta (nnz) "<< delta.nnz()<< " s(sigma) " << ssigma.nnz() << std::endl;


      /* add last row */

      M_sigma_times_delta.resize(n + 1,true); // make sure that the row has the right dimension
      M_sigma_times_delta(n) = 0;
      column(M_sigma,  n) = M_sigma_times_delta;

      /* add last column */
      s_sigma.resize(n + 1,true);
      s_sigma(n) = 0;
      row(M_sigma,  n) = s_sigma;

      //std::cout << "M(sigma) (nnz) "<<  M_sigma.nnz() <<std::endl;

      /* walk through last row to update incoming */

      const Vector& last_row = row(M_sigma,n);
      for(typename Vector::const_iterator it=last_row.begin(); it!=last_row.end(); ++it) {
	if(*it == 0.0) continue;
	a.incoming[it.index()].insert(sigma);
      }

      /* walk through last column to update outgoing */
      const Vector& last_column = column(M_sigma,n);

      for(typename Vector::const_iterator it=last_column.begin(); it!=last_column.end(); ++it) {
	if(*it == 0.0) continue;
	a.enabled[it.index()].insert(sigma);
      }
    }

    /* increase dimension of initial and accepting state vector */
    delta.resize(n+1,true);
    delta(n) = 0;


    a.setNrOfStates(n+1);
    a.initial(n) = prec_inner_prod(a.initial,delta);

    n = n + 1;
  }

  // fold into compressed form (without zero rows)
  for(TransitionMatrixIterator it=a.trans_matrix.begin();
	it != a.trans_matrix.end(); ++it) {
    int sigma (it->first);
	Matrix M_sigma ( it -> second );
	fold(M_sigma, it->second, a.trans_nz[sigma]);// minimal non-zero of a matrix
  }

  std::cout << "Backward inflate: nr of non-zeros (after)" << a.getNrOfNonzeros() << std::endl;
}

#endif
