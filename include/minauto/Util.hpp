#ifndef UTIL
#define UTIL

#include <unordered_map>
#include <unordered_set>
#include <boost/foreach.hpp>
#include "boost/shared_ptr.hpp"

#include "boost/numeric/interval.hpp"
#include "boost/numeric/interval/io.hpp"

#define foreach BOOST_FOREACH

#include<vector>
#include<list>
#include<string>
#include<set>
#include<fstream>
#include<cassert>
#include<cstring>
#include<cmath>
#include <limits>

//#include <boost/bind.hpp>

//#include <boost/rational.hpp>

#include "mpreal.h"
#include "mpreal_interval.hpp"

/*
#include <boost/thread.hpp>
#include <boost/thread/future.hpp>
*/


#include <boost/function.hpp>

#include<iostream>



namespace global_setting {
  extern bool html;
  extern int verbosity;
  extern bool interval;
  extern bool check_minimisation;
  extern bool prune;
}



template<class CharType, class CharTraits, class T, class Policies>
std::basic_ostream<CharType, CharTraits> &operator<<
  (std::basic_ostream<CharType, CharTraits> &stream,
   boost::numeric::interval<T, Policies> const &value) {
	boost::numeric::operator<< <CharType,CharTraits,T, Policies>(stream, value);
}


std::string toString(float f);
std::string toString(double f);
std::string toString(long double f);
std::string toString(int i);
std::string toString(unsigned int i);
std::string toString(const mpfr::mpreal& f);




template <typename NumericType> 
struct Precision {
	typedef typename std::pair<NumericType,NumericType> NumericTypePair;
	

	Precision(NumericType __ap, NumericType __rp = 1E-6) {
		setParameters(__ap, __rp);
	}
	
	void setParameters(NumericType __ap, NumericType __rp) {
		absolute_precision = __ap;
		relative_precision = __rp;
		rel_gap_set = std::pair<bool,bool> (false,false);
		abs_gap_set = std::pair<bool,bool> (false,false);
		tolerance = 1 ; // exp( log( absolute_precision) * 0.1);
	}

	NumericType absolute_precision;
	NumericType relative_precision;
	NumericType tolerance;
	
	bool isEqual(const NumericType& a, const NumericType& b) {
			NumericType diff(fabs(a-b));
		    if (checkThresholdAbs(diff))
		        return true;
		    NumericType relativeError = (fabs(b) > fabs(a)) ? diff/b : diff/a;
		    return checkThresholdRel(relativeError);
	}
	
	bool isZero(const NumericType& a) {
			NumericType b = fabs(a);
			return b == 0.0 || checkThresholdAbs(b);
    }
	
	bool isRelativeZero(const NumericType& a) {
			NumericType b = fabs(a);
			return b == 0.0 || checkThresholdRel(b);
	}

	
	Precision(const Precision& p) {
		*this = p;
	}

	Precision& operator=(const Precision& p) {
		absolute_precision = p.absolute_precision;
		relative_precision = p.relative_precision;
		tolerance = p.tolerance;
		rel_gap = p.rel_gap;
		abs_gap = p.abs_gap;
		rel_gap_set = p.rel_gap_set;
		abs_gap_set = p.abs_gap_set;
		return *this;
	}
	
	void printStatistics() const {
		std::cout << toString () << std::endl;
	}

	std::string toString() const {
		std::string result;

		// give diagnostics of how "close" the decisions are

		bool alc ( fabs(abs_gap.first) > 0 && 1 >= fabs(abs_gap.first * tolerance - absolute_precision) / abs_gap.first);
		bool auc ( fabs(abs_gap.second) > 0 && 1 >= fabs(abs_gap.second * tolerance - absolute_precision) / absolute_precision);

		bool rlc ( fabs(rel_gap.first) > 0 && 10 >= fabs(rel_gap.first - relative_precision) / rel_gap.first );
		bool ruc ( fabs(rel_gap.second) > 0 && 10 >= fabs(rel_gap.second - relative_precision) / relative_precision);

		std::string warning;

		if(alc)
			warning += "\nTight decision: possible false zero";
		if(auc)
			warning += "\nTight decision: possible false non-zero";

		if(rlc)
			warning += "\nTight decision: possible false relative zero";
		if(ruc)
			warning += "\nTight decision: possible false relative non-zero";

		result += "Precision: abs. precision: " 
		       + ::toString(abs_gap.first) + " <= " + ::toString (absolute_precision / tolerance) + " <= " + ::toString ( abs_gap.second )
		       + " rel. precision: " + ::toString( rel_gap.first ) +  " <= " + ::toString( relative_precision ) + " <= " + ::toString( rel_gap.second ) 
		       + warning;
		return result;
	}

public:
	bool checkThresholdAbs(const NumericType& value) {
		bool result = ( (tolerance * value) <= absolute_precision );
		if(result) {
			abs_gap.first = abs_gap_set.first ? std::max(value,abs_gap.first) : value;
			abs_gap_set.first = true;
		} else {
			abs_gap.second = abs_gap_set.second ? std::min(value,abs_gap.second) : value;
			abs_gap_set.second = true;
		}
		return result;
	}

	bool checkThresholdRel(NumericType value) {
			bool result = value <= relative_precision;
			if(result) {
				rel_gap.first = rel_gap_set.first ? std::max(value,rel_gap.first) : value;
				rel_gap_set.first = true;
			} else {
				rel_gap.second = rel_gap_set.second ? std::min(value,rel_gap.second) : value;
				rel_gap_set.second = true;
			}
			return result;
		}

	NumericTypePair rel_gap;
	NumericTypePair abs_gap;
	std::pair<bool,bool> rel_gap_set;
	std::pair<bool,bool> abs_gap_set;
};

template<class Vector>
void cleanup(Vector& v, double absolute_precision,double relative_precision) {
	typedef typename Vector::iterator VectorConstIterator;
	long double average = 0;
	unsigned count = 0;
	for(VectorConstIterator it= v.begin(); it!=v.end();++it) {
	  ++count;
	  average += fabs(*it);
	}

	average /= count; 

	for(VectorConstIterator it= v.begin(); it!=v.end();++it) {
	  if(isEqual(*it + average,average,absolute_precision,relative_precision))
	    *it = 0;
	}
}


template<class Vector>
bool isEqualVector(const Vector& v1, const Vector& v2, double absolute_precision,double relative_precision) {
	typedef typename Vector::const_iterator VectorConstIterator;
	bool z = true;
	for(VectorConstIterator it= v1.begin(); z && it!=v1.end();++it) {
	  z &= isEqual(*it,v2(it.index()),absolute_precision,relative_precision);
	}
	for(VectorConstIterator it= v2.begin(); z && it!=v2.end();++it) {
	  z &= isEqual(*it,v1(it.index()),absolute_precision,relative_precision);
	}

	return z;
}


template<class Vector>
bool vanishVector(const Vector& v1, const Vector& v2, double absolute_precision,double relative_precision) {
	typedef typename Vector::const_iterator VectorConstIterator;
	bool z = true;

/*
	std:: cout << " ==== " <<std::endl;
	std::cout << "v1 = "<< v1 << std::endl;
	std::cout << "v2 = "<< v2 << std::endl;
*/
	for(VectorConstIterator it= v1.begin(); z && it!=v1.end();++it) {
      double den (  fabs(v2(it.index())) );
	  z &= fabs(*it) / den < relative_precision;

	  if(!z) std::cout << " no: " << it.index() << " " << fabs(*it) << " " << fabs(v2(it.index())) <<  " " << relative_precision 
                       <<  " " << fabs(*it) / den << std::endl;

	}
/*
	std:: cout<< (z ? "tiny" : "significant") <<std::endl;
	std:: cout << " ==== " <<std::endl;
*/
	return z;
}


double getRandomDouble(unsigned long long K=RAND_MAX);

template<typename NumericType>
bool isEqual(NumericType& a, NumericType& b, NumericType& absolute_precision, NumericType& relative_precision) {
	NumericType diff(a-b);
    if (fabs(diff) < absolute_precision)
        return true;
    NumericType relativeError = (fabs(b) > fabs(a)) ? fabs(diff/b) : fabs(diff/a);
    return relativeError <= relative_precision;
}

template<typename Scalar>
class Comparator {
public:
	Comparator(Scalar __absolute_precision, 
			   Scalar __relative_precision)
				: abs (__absolute_precision), 
			      rel (__relative_precision), 
			      rel_gap(std::numeric_limits<Scalar>::min(),std::numeric_limits<Scalar>::max()),
			      abs_gap (std::numeric_limits<Scalar>::min(),std::numeric_limits<Scalar>::max()) {}
	
	bool operator()(Scalar a, Scalar b) {
		return check_threshold(a-b,abs,abs_gap) || check_threshold( fabs(b) > fabs(a) ? a-b/b : a-b/a, rel, rel_gap);
	}
	
	bool operator() (Scalar a) {
		return check_threshold(fabs(a),0,abs, abs_gap);
	}
	
	typedef typename std::pair<Scalar, Scalar> ScalarPair;
	
	Scalar getAbsolutePrecision() const { return abs; }
	Scalar getRelativePrecision() const { return rel; }
	ScalarPair getAbsoluteGap() const { return abs_gap; }
	ScalarPair getRelativeGap() const { return rel_gap; }

	void printStatistics() const {
		std::cout << "Comparator "<< std::endl;
		std::cout << "Absolute: " << abs_gap.first << " <= " << abs << " <= " << abs_gap.second << std::endl;
		std::cout << "Relative: " << rel_gap.first << " <= " << rel << " <= " << rel_gap.second << std::endl;
	}
	
private:
	bool check_threshold(Scalar value, Scalar threshold, ScalarPair& gap) {
		bool result = fabs(value) < threshold;
		if(result)
			gap.first = std::max(value, gap.first);
		else
			gap.second = std::min(value,gap.second);
		return result;
	}
	
	Scalar abs;
	Scalar rel;
	
	ScalarPair rel_gap;
	ScalarPair abs_gap;
};


extern void progressBar(int j, int n);

#endif
