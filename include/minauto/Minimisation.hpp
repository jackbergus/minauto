#ifndef MINIMISATION
#define MINIMISATION

#include <BLAS.hpp>

struct MinimisationSettings {

  MinimisationSettings () : method (arnoldi), direction(both), prec(false), dense(false) {}

  enum Method {
    arnoldi,
    householder,
    compare
  } method;

  enum Direction {
	  forward,
	  backward,
	  both
  } direction;

  bool prec;   // use long doubles
  bool dense;  // use dense matrices
  bool lz;
  bool normalise;

  bool pivot;
  bool reortho;
  bool classical; // classical Gram-Schmidt
};



template <typename Automaton>
struct Witness {
  typedef typename Automaton::Matrix Matrix;
  typedef typename Automaton::Weight Scalar;

  Witness () : cut (false), orthogonality(0) {}

  /* basis of the row or column space */
  Matrix F;

  /* language tree */
  WordTree tree;

  /* statistics */
  long iterations;
  long time;                           // time in seconds

  bool cut;
  std::vector<std::pair<int,int> > cutpoints ; // cut points in the language tree

  Scalar orthogonality;
  std::unordered_map<int,Scalar> commutativity; // commutativity of reduction

  // print a report
  std::string toString() const {
	std::string result;

	result += "Reduction from " + (::toString((unsigned int) F.size2())) + " to " + (::toString((unsigned int)F.size1())) + " states \n";
	result += "orthogonality "+ ::toString(orthogonality) +  " \n";

	if( cut ) {
		std::cout << " === Warning: the result might be unsafe. === " << std::endl;
	}

	return result;
  }

};


template<typename Vector>
bool checkInclusion(const std::vector<Vector>& V, const Vector& w) {

	std::string result;

	std::cout << "Linear Programming problem" << std::endl;

	typedef typename Vector::const_iterator VectorConstIterator;

	/*
	 * generates an LP with unknowns X={x0, ..., xn }
	 * where
	 * - V = {V^(0), ... , V^(n)}
	 * - n = size(base) - 1
	 * - m = dim V^(i) - 1
	 *
	 * Minimise
	 *   obj:
	 * Subject To
	 *   V^(0)_0 * x0 + ... + V^(n)_1 * xn <= w_0
	 *                ...
	 *   V^(0)_m * x0 + ... + V^(n)_m * xn <= w_m
	 * Bounds
	 *   x0 >= 0
	 *    ...
	 *   xn >= 0
	 * End
	 */
	result += "Minimise\n";
	result += "  obj:\n";
	result += "Subject To\n";

	if(V.size() > 0) {

		unsigned m = V[0].size();

		std::string constraints[m];

		for(unsigned i=0; i<V.size();++i) {
			for(VectorConstIterator it=V[i].begin();it != V[i].end();++it) {
				unsigned j = it.index();
				constraints[j] += " + " + toString(*it) + " * x" + toString(i);
			}
		}
		for(unsigned i=0; i<m;++i) {
			if(constraints[i].size()>0) {
				result += constraints[i] + " >= " + toString(w(i))+ "\n";
			} else if (w(i)!= 0){
				result += "  0 >= " + toString(w(i))+ "\n";
			}
		}
	}
	result += "Bounds\n";
	for(unsigned i=0; i<V.size();++i) {
		result += "  x" + toString(i) + " >= 0\n";
	}

	result += "End\n";

	std::cout << result << std::endl;

	return true;
}



template <typename Automaton>
void checkLargerZero(Automaton& a, Precision<typename Automaton::Weight>& prec) {
  typedef typename Automaton::Weight Scalar;
  typedef typename Automaton::Vector Vector;
  typedef typename Automaton::Matrix Matrix;



  unsigned n = a.getNrOfStates();

  std::vector<Vector> history;

  // compute the vector space span(F) = { alpha M(w) | word w}


  std::unordered_map<int, std::vector<std::pair<unsigned, Vector> > > combinations;

  Vector init(a.initial);

  history.push_back(init);

  for (unsigned j = 1, l = 0; l < j; ++l) {
	  assert(l<history.size());
    const Vector& current = history[l];

    std::unordered_set<int> enabled;
    a.get_enabled(current,enabled);

    foreach(int sigma, enabled) {
      Vector next (n);
      std::cout << "Image under of " << current << " under " << alphabet[sigma] << std::endl;
      a(sigma, next, current, true);


      if(inner_prod(next,a.accepting) < 0) {
    	  std::cout << "Not smaller zero " << std::endl;
    	  return;
      }

      std::cout << "next " << next << std::endl;
      std::cout << "accepting " << a.accepting << std::endl;
      std::cout << "next' * accepting " << inner_prod(next,a.accepting) << std::endl;

      Scalar norm(norm_2(next));
      if(norm > prec.absolute_precision ) {
    	  checkInclusion(history,next);

          std::cout << "Storing result in history " << std::endl;

    	  history.push_back(next);

    	  ++j;
      }
      std::cout << " === " << std::endl;
    }
  }




}

template <typename Automaton>
bool isContained(typename Automaton::Vector& v, typename Automaton::IntervalVector& iv) {
	typedef typename Automaton::Interval Interval;
	typedef typename Automaton::Vector Vector;
	bool result = true;
	
	for(typename Vector::iterator it = v.begin(); it!=v.end(); ++it) {
		Interval& inter(iv(it.index()));
		result &= (*it >= inter.lower()) && (*it <= inter.upper());		
	}
	return result;
	
}







template <typename Automaton, bool forward>
void buildTransitions(	Automaton& a,
			Automaton& b,
			std::unordered_map<int, std::vector<std::pair<unsigned, typename Automaton::Vector> > >& combinations,
			Precision<typename Automaton::Weight>& prec,
			Witness<Automaton>& wit)
{
	typedef typename Automaton::Weight Scalar;
	typedef typename Automaton::Vector Vector;
	typedef typename Vector::const_iterator VectorConstIterator;
	typedef typename Automaton::Matrix Matrix;

	typedef typename Automaton::TransitionMatrix TransitionMatrix;
	typedef typename TransitionMatrix::iterator TransitionMatrixIterator;


	const Matrix& F = wit.F;
	unsigned nf = F.size1();


	if(global_setting::prune) {
		filter(b.initial,prec);
		filter(b.accepting,prec);
	}


	Scalar max_angle(0);

	int orth_i = 0, orth_j = 0;


  if(global_setting::verbosity > 1) {
	  std::cout << "Orthogonality check " ;
	  for(unsigned i=0; i<F.size1(); ++i) {
		  for(unsigned j=i+1; j<F.size1();++j) {
			  Scalar angle(fabs(inner_prod(row(F,i),row(F,j))));
			  if( angle > max_angle) {
				  max_angle = angle;
				  orth_i = i;
				  orth_j = j;
			  }
		  }
	  }
	  wit.orthogonality = std::max(wit.orthogonality,max_angle);

	  std::cout << " largest covariance " << orth_i << " " << orth_j << std::endl;
  }

  /*
	if(global_setting::check_minimisation) {
		Vector initial(prec_prod(b.initial,F) - a.initial);
		Scalar initial_diff = norm_inf(initial);
		int arg_initial_diff = index_norm_inf(initial);
		Vector accepting (prec_prod(trans(F),b.accepting) - a.accepting);
		Scalar accepting_diff = norm_inf(accepting);
		int arg_accepting_diff = index_norm_inf(accepting);
		std::cout << "initial difference " << initial_diff << " at " << arg_initial_diff
				  <<" accepting difference " << accepting_diff << " at " << arg_accepting_diff << std::endl;
	}
	*/

	// ... enumerate the different letters and go through
	for(TransitionMatrixIterator it=a.trans_matrix.begin();
		it != a.trans_matrix.end(); ++it) {
	  int sigma (it->first);
	  Matrix Mf_sigma(nf,nf);
	  Mf_sigma.clear();

	  std::vector< std::pair<unsigned,Vector> >& combinations_sigma(combinations[sigma]);

	  if(forward) {
		  for(unsigned ri=0; ri!=combinations_sigma.size(); ++ri) {
			unsigned source = combinations_sigma[ri].first;
			Vector& row = combinations_sigma [ri].second;

			row = project(row,boost::numeric::ublas::range(0,nf));

			if(global_setting::prune)
				filter(row,prec);

			b.enabled[source].insert(sigma);

			for(VectorConstIterator rit = row.begin(); rit!=row.end(); ++rit) {
			  if(*rit == 0) continue;

			  unsigned target = rit.index();
			  Mf_sigma(source,target) = *rit;
			  b.incoming[target].insert(sigma);
			}
		  }
	  } else {
		 for(unsigned ri=0; ri!=combinations_sigma.size(); ++ri) {
	        	unsigned target = combinations_sigma[ri].first;
		        Vector& row = combinations_sigma [ri].second;

			row = project(row,boost::numeric::ublas::range(0,nf));

			if(global_setting::prune)
				filter(row,prec);

			b.incoming[target].insert(sigma);

		        for(VectorConstIterator rit = row.begin(); rit!=row.end(); ++rit) {
        			if(*rit == 0) continue;

			          unsigned source = rit.index();
				  Mf_sigma(source,target) = *rit;
				  b.enabled[source].insert(sigma);
			}
	  	}
	  }
	  
    /*
	  if(global_setting::check_minimisation) {
		  Matrix M;
		  unfold(it->second, a.trans_nz[sigma], M);

		  Scalar diff = 0;
		  Scalar fdiff = 0;
		  int i = 0, j = 0;

		  if(forward) {
			  Matrix Diff (prod(Mf_sigma,F) - prod(F,M));
			  fdiff = norm_inf(Diff);

			  Scalar max = 0;
			  for(typename Matrix::const_iterator1 it1 = Diff.begin1(); it1!=Diff.end1(); ++it1) {
				  for(typename Matrix::const_iterator2 it2 = it1.begin() ; it2 != it1.end(); ++it2) {
					  Scalar d(fabs(*it2));
					  max = std::max(max,d);
				  }
			  }

			  std::cout << "max " << max << std::endl;

			  for(unsigned r=0; r<Diff.size1(); ++r) {

				  Scalar rdiff (norm_inf(row(Diff,r)));
				  if(diff < rdiff) {
					  diff = rdiff;
					  i = r;
					  j = index_norm_inf(row(Diff,r));
				  }
			  }
		  } else {
			  Matrix Diff (prod(trans(F),Mf_sigma) - prod(M,trans(F)));

			  fdiff = norm_inf(Diff);
			  for(unsigned r=0; r<Diff.size1(); ++r) {

				  Scalar rdiff (norm_inf(row(Diff,r)));
				  if(diff < rdiff) {
					  diff = rdiff;
					  i = r;
					  j = index_norm_inf(row(Diff,r));
				  }
			  }
		  }

		  wit.commutativity[ sigma ] = diff;
		  std::cout << "difference " << alphabet[sigma] << " max_i,j |M' F - F M| " <<  fdiff << " max entry " << diff << " at (" <<i<< "," << j <<")" << std::endl;
	  }
	  */

	  fold(Mf_sigma,b.trans_matrix[sigma],b.trans_nz[sigma]);
	}
	//b.reachabilityAnalysis(!forward);
}



/* Compute the reflection vector
 * pre:  x contains a non-zero vector
 * post: x contains the reflection vector and scalar its norm
 *
 * See: http://z.cs.utexas.edu/wiki/LA.wiki/Householder_Transformations/Unblocked
 */
template <class Vector, class Scalar>
void getReflector(Vector& x, Scalar& scalar) {
	boost::numeric::ublas::vector_range<Vector> x2(x, boost::numeric::ublas::range(1,x.size()));

	/*
	x(0) += sign(x(0)) * norm_2(x);
	x /= norm_2(x);
        scalar = 0.5;
	*/
	/*

	Scalar chi2 = prec_inner_prod(x2,x2);

	Scalar pivot(x(0) + x(0) >= 0 ? norm_2(x) : -norm_2(x) );
	x(0) = 1;
	x2 /= pivot;
	chi2 /= pivot * pivot;
	scalar = (1 + chi2) * 0.5;
	*/
	
	x(0) += (x(0) >=0 ? 1 : -1) * norm_2(x);
	scalar = prec_inner_prod(x,x) * 0.5;
}

template <typename Vector>
void reflectVector(const Vector& reflector,
				   const typename Vector::value_type& scalar,
				   Vector& v,
				   unsigned j = 0) {
	boost::numeric::ublas::vector_range<Vector> vrange(v, boost::numeric::ublas::range(j,v.size()));
	typename Vector::value_type inner(prec_inner_prod(vrange,reflector)	);

	vrange -= (inner/scalar) * reflector;

	/*vrange.minus_assign(prec_prod(outer_prod(reflector,reflector), vrange) / scalar);*/
}




template <typename Vector, bool forward>
void reflectVector(const std::vector<Vector>& reflectors,
				   const Vector& scalars,
		           Vector& v) {
	if(forward) {
		for(unsigned j=0; j<reflectors.size(); ++j) {
			reflectVector(reflectors[j], scalars(j), v, j);
		}
	} else {
		for(int j=reflectors.size() -1; j>=0 ; --j) {
			reflectVector(reflectors[j], scalars(j), v, j);
		}
	}
}


template <typename Automaton, bool forward>
void Householder(Automaton& a,
		 Automaton& b,
		 MinimisationSettings& s,
		 Precision<typename Automaton::Weight>& prec,
		 Witness<Automaton>& wit) {
  typedef typename Automaton::Weight Scalar;
  typedef typename Automaton::Interval Interval;
  typedef typename Automaton::Vector Vector;
  typedef typename Automaton::ZeroVector ZeroVector;
  typedef typename Automaton::Matrix Matrix;
  typedef typename Automaton::IdentityMatrix IdentityMatrix;
  typedef typename boost::numeric::ublas::vector_range< Vector > VectorRange;
  typedef typename boost::numeric::ublas::unit_vector<Scalar> UnitVector;
  typedef typename Vector::const_iterator VectorConstIterator;
  typedef typename Vector::iterator VectorIterator;

  unsigned n = a.getNrOfStates();

  // compute the vector space span(F) = { alpha M(w) | word w}
  std::vector<Vector> basis;
  std::vector<Vector> reflectors;
  Vector scalars(n, 0);

  std::unordered_map<int, std::vector<std::pair<unsigned, Vector> > > combinations;

  Vector start(forward ? a.initial : a.accepting);

  start /= norm_2(start);

  Vector reflector_start (start);
  Scalar scalar(0);

  getReflector<Vector,Scalar>(reflector_start,scalar);
  scalars(0) = scalar;

  Vector start_basis(n, 0);
  start_basis(0) = 1;

  reflectVector<Vector>(reflector_start,scalar,start_basis);

  Scalar sign( prec_inner_prod(start_basis,start) );

  reflectors.push_back(reflector_start);

  basis.push_back(start_basis);

  unsigned j = 1;

  for (unsigned l = 0; l < j; ++l) {
	  std::unordered_set<int> letters;
	  a.getLetters(basis[l],letters,forward);

	  progressBar(j,n);
	  foreach(int sigma, letters) {

      Vector next (n, 0);

	    a(sigma, next, basis[l], forward);


	    VectorRange next_range(next,boost::numeric::ublas::range(j,n));
	    reflectVector<Vector,true>(reflectors,scalars,next);

	    bool stop ( prec.isZero(norm_inf(next_range) ) || prec.isRelativeZero(norm_2(next_range) / norm_2(next)));

	    if( !stop ) {
		    if(j >=n ) { 
		      wit.cut = true; // this indicates a numerical problem
		    } else {

		      Vector reflector(next_range);

		      getReflector<Vector,Scalar>(reflector,scalar);
          scalars(j) = scalar;
		      reflectVector<Vector>(reflector,scalars(j),next,j);

		      reflectors.push_back(reflector);


		      /* new basis vector */

          Vector basis_vector(n, 0);
		      basis_vector(j) = 1;
		      reflectVector<Vector,false>(reflectors, scalars, basis_vector);

		      basis.push_back(basis_vector);
		      ++j;
		    }
	    } 
	    combinations[sigma].push_back(std::pair<unsigned,Vector>(l,next));
	  }
	
  }


  // compute the new transition matrices
  unsigned nf = basis.size();

  b.setNrOfStates(nf);

  // initial and accepting vector

  wit.F.resize(nf,n,false);

  for(unsigned i=0; i<basis.size(); ++i) {
	  row(wit.F,i) = basis[i];
  }


  if(forward) {
	  b.accepting = prec_prod(wit.F,a.accepting);

	  b.initial.clear();
	  b.initial.resize(nf,false);
	  b.initial(0) = sign * norm_2(a.initial);
  } else {
	  b.accepting.resize(nf,false);
	  b.accepting(0) = sign * norm_2(a.accepting);
	  b.initial = prec_prod(wit.F,a.initial);
  }

  /*
  if(forward) {
	  b.initial.resize(nf,false);
	  b.initial.clear();
	  b.initial(0) = sign * norm_2(a.initial);
	  b.accepting = a.accepting;
	  reflectVector<Vector,true>(reflectors, scalars, b.accepting);
	  b.accepting.resize(nf,true);
  } else {
	  b.accepting.resize(nf,false);
	  b.accepting.clear();
	  b.accepting(0) = sign * norm_2(a.accepting);
	  b.initial = a.initial;
	  reflectVector<Vector,true>(reflectors, scalars, b.initial);
	  b.initial.resize(nf,true);
  }
  */
  buildTransitions<Automaton,forward>(a, b, combinations, prec, wit);
}






template <typename Automaton, bool forward>
void Arnoldi(Automaton& a, 
	     Automaton& b, 
	     MinimisationSettings& s, 
	     Precision<typename Automaton::Weight>& prec,
  	     Witness<Automaton>& wit) {
  typedef typename Automaton::Weight Scalar;
  typedef typename Automaton::Vector Vector;
  typedef typename Automaton::ZeroVector ZeroVector;
  typedef typename Automaton::Matrix Matrix;

  unsigned n = a.getNrOfStates();

  GramSchmidt<Automaton> gs(n,prec,true,s.pivot,s.reortho,s.classical);

  // compute the vector space span(F) = { alpha M(w) | word w}
  


  std::unordered_map<int, std::vector<std::pair<unsigned, Vector> > > combinations;

  Vector start(forward ? a.initial : a.accepting);

  start /= norm_2(start);

  gs.addToBasis(start);

  Vector next (n);

  for (unsigned j = 1, l = 0; l < j; ++l) {
    std::unordered_set<int> letters;
    a.getLetters(gs[l],letters,forward);
    
    progressBar(j,n);
    
    if(global_setting::verbosity > 1) {
      std::cout << "norm " << norm_inf(gs[l]) << "\n";
      std::cout << "<<<<< \nBase [" << l << "] " << gs[l] << " >>>> \n";
    }

    Vector row(std::max(j+1,n));

    foreach(int sigma, letters) {

      next.clear();
      a(sigma, next, gs[l], forward);

      if(global_setting::prune)
      	      filter<Vector>(next, prec.absolute_precision);

      if(global_setting::verbosity > 1) {
      	std::cout << "====\nletter " << alphabet[sigma] << " j " << j << std::endl;

      	std::cout << "A(w) = " << next << std::endl;
      }

      row.clear();
      gs.project(next,row);
      
      if(global_setting::prune)
      {
      	      filter<Vector>(row,prec.absolute_precision);
      	      filter<Vector>(next, prec.absolute_precision);
      }

      if(global_setting::verbosity > 1) {
      	std::cout << "=> \n" << next << std::endl;
      }

      Scalar norm(norm_2(next));

      if( norm > prec.absolute_precision ) {
    	  if( j >= n ) { // this indicates a numerical problem
    	    wit.cut = true;
    	  } else {
    	  	row(j) = norm;
    	  	++j;
    	  	next /= norm;
    	   	gs.addToBasis(next);
    	  }
      }

      combinations[sigma].push_back(std::pair<unsigned,Vector>(l,row));
    }
  }
  
  // compute the new transition matrices
  unsigned nf = gs.getDimension();

  b.setNrOfStates(nf);
  
  gs.printStatistics(global_setting::verbosity);

  wit.F.resize(nf,n,false);

  for(unsigned i=0; i<nf; ++i) {
	  row(wit.F,i) = gs[i];
  }

  if(forward) {
	  b.accepting = prec_prod(wit.F,a.accepting);

	  b.initial.clear();
	  b.initial.resize(nf,false);
	  b.initial(0) = norm_2(a.initial);
  } else {
	  b.accepting.resize(nf,false);
	  b.accepting(0) = norm_2(a.accepting);
	  b.initial = prec_prod(wit.F,a.initial);
  }

  buildTransitions<Automaton,forward>(a, b, combinations, prec, wit);
}

template <typename Automaton, bool forward>
void Reduce(MinimisationSettings::Method& m,
	  Automaton& a, 
	  Automaton& b, 
	  MinimisationSettings& s, 
	  Precision<typename Automaton::Weight>& prec,
  	  Witness<Automaton>& wit) {

	  typedef typename Automaton::Matrix Matrix;

	  typedef typename Automaton::TransitionMatrix::const_iterator TransitionMatrixIterator;

	  typedef typename Automaton::Vector Vector;

  switch(s.method) {
	case MinimisationSettings::arnoldi:

		std::cout << "Gram-Schmidt " << (forward ? "forward" : "backward" ) << std::endl;

		Arnoldi<Automaton,forward>(a,b,s,prec,wit);

		std::cout << wit.toString() << std::endl;
		break;
	case MinimisationSettings::householder:

		std::cout << "Householder " << (forward ? "forward" : "backward" ) << std::endl;

		Householder<Automaton,forward>(a,b,s,prec,wit);

		std::cout << wit.toString() << std::endl;
		break;
        case MinimisationSettings::compare:
	  {
	    Witness<Automaton> wit_arno;
	    Automaton b_arno;
	    Arnoldi<Automaton,forward>(a,b_arno,s,prec,wit_arno);
	    Householder<Automaton,forward>(a,b,s,prec,wit);

	    typedef typename Automaton::Matrix Matrix;

	    /* compute the sign-flip matrix
	     * Flip is a square diagonal matrix with unit entries unless numerical errors occur
	     */
	    Matrix Flip = prod( wit_arno.F , boost::numeric::ublas::trans(wit.F));
	    
	    /* check if the automata are equal up to signs
	     * alpha' == alpha * Flip ?
	     * eta'   == Flip * eta ?
	     * foreach(sigma in alphabet)
	     *   M'(sigma) == Flip * M(sigma) * Flip ?
	     */
	    
	    Vector initial_diff(b.initial - prod(b_arno.initial, wit.F));
	    Vector accepting_diff(b.accepting - prod(boost::numeric::ublas::trans(Flip),b_arno.accepting));

		for(TransitionMatrixIterator it=a.trans_matrix.begin();
			it != a.trans_matrix.end(); ++it) {
		  int sigma (it->first);

		  /*
		   * 1. unfold both transition matrices
		   * 2. compare them
		   */

		  Matrix M_house, M_arno;

		  unfold(M_arno, b_arno.trans_nz[sigma], b_arno.trans_matrix[sigma]);
		  unfold(M_house, b.trans_nz[sigma], b.trans_matrix[sigma]);

		  // compute the difference

	    }


	  }
	  
	        break;
	default:
		std::cout << "Minimisation method not supported" << std::endl;
		break;
  }
}


template <typename Automaton>
int minimise(Automaton& a, Automaton& b, MinimisationSettings& s, Precision<typename Automaton::Weight>& prec) {

  typedef typename Automaton::Matrix Matrix;
  Timer min_timer;

  Witness<Automaton> wit1;
  Witness<Automaton> wit2;

  min_timer.Start();

  switch(s.direction) {
  	  case MinimisationSettings::both:
	  {
		  Automaton c;
		  Reduce<Automaton,true>(s.method,a,c,s,prec,wit1);
		  Reduce<Automaton,false>(s.method,c,b,s,prec,wit2);
	  }
  		  break;
  	  case MinimisationSettings::forward:
  		  Reduce<Automaton,true>(s.method,a,b,s,prec,wit1);
  		  break;
	  case MinimisationSettings::backward:
		  Reduce<Automaton,false>(s.method,a,b,s,prec,wit1);
	  break;
  }

  prec.printStatistics();

  if(s.lz) {
	  checkLargerZero(b, prec);
  }

  min_timer.Stop();
  return 0;
}

#endif
