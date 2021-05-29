#ifndef GRAMSCHMIDT
#define GRAMSCHMIDT

template<typename Automaton>
class GramSchmidt {
public:
	typedef typename Automaton::Matrix Matrix;
	typedef typename Automaton::Vector Vector;
	typedef typename Automaton::Weight Scalar;
	typedef typename Automaton::Interval Interval;
	typedef typename Automaton::IntervalVector IntervalVector;

	GramSchmidt(unsigned dim, 
		    Precision<typename Vector::value_type>& __prec, 
		    bool n = false,
		    bool __pivot = false,
        bool __reortho = false,
		    bool __classical = false) :
		prec(__prec), 
		dim(0), 
		nzrows(dim), 
		cutoffs(0), 
		normalise(n),
		pivot(__pivot),
		reortho (__reortho),
                classical (__classical) {
	}

	void project(Vector& u) {
		Scalar norm_u(norm_2(u));
		u /= norm_2(u);

		timer.Start();

		getRelevantBaseVectors(u, axes);

		
		// classical projection
		if(classical)
		{
		  Vector sum(u.size(), 0);
		  foreach(int i, axes) {
		    const Vector& b(basis[i]);
		    Scalar scalar = normalise ? prec_inner_prod(b,u) : prec_inner_prod(b,u) / inner_prods[i];
		    sum.plus_assign(scalar * b);
		  }
		  u.minus_assign(sum);
		}
		else
		{
		  // more stable projection
		  foreach(int i, axes) {
		    const Vector& b(basis[i]);
		    Scalar scalar = normalise ? prec_inner_prod(b,u) : prec_inner_prod(b,u) / inner_prods[i];
		    u.minus_assign(scalar * b);
		  }
		}

		//filterVector(u, prec.relative_precision);

		if ( prec.isZero(norm_inf(u)) || prec.isRelativeZero(norm_2(u) / norm_u) ) {
			u.clear();
		}

		timer.Stop();
	}

	void project(Vector& u, Vector& scalars) {
		timer.Start();

		Scalar norm_u(norm_2(u));
		
		getRelevantBaseVectors(u, axes);

		/* Kahan summation */

		/*
		Vector c(u.size());
		foreach(int i, axes) {
			const Vector& b(basis[i]);
			Scalar scalar = normalise ? prec_inner_prod(b,u) : (prec_inner_prod(b,u) / inner_prods[i]);
			Vector t(u - (scalar * b + c));
			c = (t - u ) + (scalar * b + c);
			scalars(i) = scalar;
			u = t;
		}
		*/

		// classical projection
		if(classical)
		{
		  Vector sum(u.size(),0);
		  foreach(int i, axes) {
		    const Vector& b(basis[i]);
		    Scalar scalar = normalise ? prec_inner_prod(b,u) : (prec_inner_prod(b,u) / inner_prods[i]);
		    scalars(i) = scalar;
		    sum.plus_assign(scalar * b);
		  }
		  u.minus_assign(sum);
		}
		else
		{
		  foreach(int i, axes) {
		    const Vector& b(basis[i]);
		    Scalar scalar = normalise ? prec_inner_prod(b,u) : (prec_inner_prod(b,u) / inner_prods[i]);
		    scalars(i) = scalar;
		    u -= scalar * b;
		  }
		}
		if(reortho) {
			foreach(int i, axes) {
				const Vector& b(basis[i]);
				Scalar scalar = normalise ? prec_inner_prod(b,u) : (prec_inner_prod(b,u) / inner_prods[i]);
				scalars(i) += scalar;
				u -= scalar * b;
			}
		}

		//std::cout << "GramSchmidt summation error " << norm_inf(c) << std::endl;
		if (  prec.isZero(norm_inf(u)) ||  prec.isRelativeZero(norm_2(u) / norm_u) ) {
			u.clear();
		}

		timer.Stop();
	}

	bool addToBasis(Vector& q) {

		bool independent = true; // is q linearly independent from basis?

		Scalar norm = norm_inf(q);

		if (norm == 0) {
			independent = false;
		} else {
			timer.Start();
			
			if (dim >= q.size()) {
				++cutoffs;
			}
			
			basis.push_back(q);
	
			if (!normalise)
				inner_prods.push_back(prec_inner_prod(q, q));

			updateNzrows(q);
			++dim;

			timer.Stop();
		}
		return independent;
	}

	/** Statistics */

	long getMem() const {
		long result = 0;
		for(unsigned i=0; i<dim; ++i)
			result += mem(basis[i]);
		return result + nzrows.size() * sizeof(int);
	}

	unsigned getDimension() const {
		return dim;
	}

	unsigned getCutOffs() const {
		return cutoffs;
	}

	double getTime() const {
		return timer.Read();
	}

	const Vector& operator[] (unsigned i) const {
		return basis[i];
	}

	void printStatistics(int verbosity = 0) const {
		std::cout << (classical ? "Classical "  : " ") << "Gram Schmidt :" << std::endl;
		std::cout << "   time " << getTime() << std::endl;
	}

private:
	typedef typename Matrix::iterator1 MatrixIterator1;
	typedef typename Matrix::iterator2 MatrixIterator2;
	typedef typename Matrix::const_iterator1 MatrixConstIterator1;
	typedef typename Matrix::const_iterator2 MatrixConstIterator2;
	typedef typename Vector::const_iterator VectorConstIterator;
	typedef typename Vector::iterator VectorIterator;
	typedef std::unordered_set<int> Set;

	std::vector<Vector> basis;
	Precision<typename Vector::value_type>& prec;
	
	std::vector<Scalar> inner_prods;
	unsigned dim;
	std::vector<std::vector<int> > nzrows;

	// diagnostics
	unsigned cutoffs;


	// settings
	bool normalise;
  bool pivot;            // use pivoting
	bool reortho;          // use re-orthogonalisation
  bool classical;

	std::vector<int> axes;

	Timer timer;

	void updateNzrows(const Vector& u) {
		// iterate over the non-zeros of the vector
		//std::cout<<"updateNzrows "<< u << std::endl;
		for (VectorConstIterator it = u.begin(); it != u.end(); ++it) {
			if (*it == 0.0)
				continue;
			unsigned index = it.index();
			nzrows[index].push_back(dim);
		}
	}

	void getNzrows(const Vector& u, Set& rows) {
		for (VectorConstIterator it = u.begin(); it != u.end(); ++it) {
		  if (*it == 0.0)
		    continue;
		  unsigned index = it.index();
		  rows.insert(nzrows[index].begin(), nzrows[index].end());
		}
	}

	void getRelevantBaseVectors(const Vector& u, std::vector<int>& axes) {

		Set rows;
		Scalar scalar;

		axes.clear();

		// determine which base vectors share positions with non-zero entries
		getNzrows(u, rows);

		axes.reserve(rows.size());
		std::unordered_map<int, Scalar> scalars;

		// compute scalar product along all axes
		for (Set::const_iterator it = rows.begin(); it != rows.end(); ++it) {			
			const Vector& b(basis[*it]);
			scalar = fabs(normalise ? prec_inner_prod(b, u) : (prec_inner_prod(b, u) / inner_prods[*it]));
			scalars[*it] =  scalar;

			axes.push_back(*it);
		}

		// sort by ascending absolute value
	
		if(pivot) {
			compareByScalar comparator(scalars);
			std::sort(axes.begin(), axes.end(), comparator);
		}
	}


	class compareByScalar {
	public:

		std::unordered_map<int, Scalar>& scalars;

		compareByScalar(std::unordered_map<int, Scalar>& __scalars) :
			scalars(__scalars) {
		}

		bool operator()(int i1, int i2) {
			return scalars[i1] < scalars[i2];
		}
	};

};

#endif
