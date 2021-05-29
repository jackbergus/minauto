/*
 * mpreal_interval.hpp
 *
 *  Created on: Feb 23, 2012
 *      Author: Dr. Bjorn Wachter, University of Oxford
 */

#ifndef MPREAL_INTERVAL_HPP_
#define MPREAL_INTERVAL_HPP_



namespace boost {
  namespace numeric {
    namespace interval_lib {

	/* Rounding requirements */

	template<> class rounded_math<mpfr::mpreal> {
	public:
	  // default constructor, destructor
	  rounded_math<mpfr::mpreal>() {}
	  ~rounded_math<mpfr::mpreal>() {}

	  const static mp_rnd_t up = MPFR_RNDU;
	  const static mp_rnd_t down = MPFR_RNDD;

	  typedef mpfr::mpreal T;


	  // mathematical operations
	  static T add_down(const T& a, const T& b) { return add(a,b,MPFR_RNDD); }
	  static T add_up  (const T& a, const T& b) { return add(a,b,MPFR_RNDU); }
	  static T sub_down(const T& a, const T& b) { return sub(a,b,MPFR_RNDD); }
	  static T sub_up  (const T& a, const T& b) { return sub(a,b,MPFR_RNDU); }
	  static T mul_down(const T& a, const T& b) { return mul(a,b,MPFR_RNDD); }
	  static T mul_up  (const T& a, const T& b) { return mul(a,b,MPFR_RNDU); }
	  static T div_down(const T& a, const T& b) { return div(a,b,MPFR_RNDD); }
	  static T div_up  (const T& a, const T& b) { return div(a,b,MPFR_RNDU); }
	  static T sqrt_down(const T& a)     { return sqrt(a,MPFR_RNDD); }
	  static T sqrt_up  (const T& a)     { return sqrt(a,MPFR_RNDU); }
	  static T exp_down(const T& a)      { return exp(a,MPFR_RNDD); }
	  static T exp_up  (const T& a)      { return exp(a,MPFR_RNDU); }
	  static T log_down(const T& a)      { return log(a,MPFR_RNDD); }
	  static T log_up  (const T& a)      { return log(a,MPFR_RNDU); }
	  static T cos_down(const T& a)      { return cos(a,MPFR_RNDD); }
	  static T cos_up  (const T& a)      { return cos(a,MPFR_RNDU); }
	  static T tan_down(const T& a)      { return tan(a,MPFR_RNDD); }
	  static T tan_up  (const T& a)      { return tan(a,MPFR_RNDU); }
	  static T asin_down(const T& a)     { return asin(a,MPFR_RNDD); }
	  static T asin_up  (const T& a)     { return asin(a,MPFR_RNDU); }
	  static T acos_down(const T& a)     { return acos(a,MPFR_RNDD); }
	  static T acos_up  (const T& a)     { return acos(a,MPFR_RNDU); }
	  static T atan_down(const T& a)	 { return atan(a,MPFR_RNDD); }
	  static T atan_up  (const T& a)     { return atan(a,MPFR_RNDU); }
	  static T sinh_down(const T& a)     { return sinh(a,MPFR_RNDD); }
	  static T sinh_up  (const T& a)     { return sinh(a,MPFR_RNDU); }
	  static T cosh_down(const T& a)     { return cosh(a,MPFR_RNDD); }
	  static T cosh_up  (const T& a)     { return cosh(a,MPFR_RNDU); }
	  static T tanh_down(const T& a)     { return tanh(a,MPFR_RNDD); }
	  static T tanh_up  (const T& a)     { return tanh(a,MPFR_RNDU); }
	  static T asinh_down(const T& a)    { return asinh(a,MPFR_RNDD); }
	  static T asinh_up  (const T& a)    { return asinh(a,MPFR_RNDU); }
	  static T acosh_down(const T& a)    { return acosh(a,MPFR_RNDD); }
	  static T acosh_up  (const T& a)    { return acosh(a,MPFR_RNDU); }
	  static T atanh_down(const T& a)    { return atanh(a,MPFR_RNDD); }
	  static T atanh_up  (const T& a)    { return atanh(a,MPFR_RNDU); }
	  static T median(const T& a, T b)   { return agm(a,b); }
	  static T int_down(const T& a)      { return rint_floor(a,MPFR_RNDD); }
	  static T int_up  (const T& a)      { return rint_ceil(a,MPFR_RNDU); }
	  // conversion functions
	  template<typename U>
	  T conv_down(const U& u)     { return T(u,mpfr::mpreal::default_prec,MPFR_RNDD); }
	  template<typename U>
	  T conv_up  (const U& u)     { return T(u,mpfr::mpreal::default_prec,MPFR_RNDU); }
	  // unprotected rounding class
	  typedef rounded_math<mpfr::mpreal> unprotected_rounding;
	};

	template<> class checking_base<mpfr::mpreal> {
	public:
		static mpfr::mpreal pos_inf() {
			return mpfr::mpreal().setInf();
		}
		static mpfr::mpreal neg_inf() {
			return mpfr::mpreal().setInf();
		}

		static mpfr::mpreal nan() {
			return mpfr::mpreal().setNan();
		}
		static bool is_nan(const mpfr::mpreal& x) {
			return isnan (x);
		}

		static mpfr::mpreal empty_lower() {
			return mpfr::mpreal().setNan();
		}

		static mpfr::mpreal empty_upper() {
			return mpfr::mpreal().setNan();
		}

		static bool is_empty(const mpfr::mpreal& x, const mpfr::mpreal& y) {
			return isnan (x) && isnan (y);
		}
	};

	namespace user {
		inline bool is_zero(const mpfr::mpreal& x) { return iszero (x); }
		inline bool is_neg (const mpfr::mpreal& x) { return sgn (x) < 0; }
		inline bool is_pos (const mpfr::mpreal& x) { return sgn (x) > 0; }
	}

	namespace constants {
		template<> inline mpfr::mpreal pi_half_lower<mpfr::mpreal>() {
			return mpfr::const_pi (mpfr::mpreal::default_prec, MPFR_RNDD)  / 2.0;
		}
		template<> inline mpfr::mpreal pi_half_upper<mpfr::mpreal>() {
			return mpfr::const_pi (mpfr::mpreal::default_prec, MPFR_RNDU)  / 2.0;
		}
		template<> inline mpfr::mpreal pi_lower<mpfr::mpreal>() {
			return mpfr::const_pi (mpfr::mpreal::default_prec, MPFR_RNDD);
		}
		template<> inline mpfr::mpreal pi_upper<mpfr::mpreal>() {
			return mpfr::const_pi (mpfr::mpreal::default_prec, MPFR_RNDU);
		}
		template<> inline mpfr::mpreal pi_twice_lower<mpfr::mpreal>() {
			return mpfr::const_pi (mpfr::mpreal::default_prec, MPFR_RNDD)  * 2.0;
		}
		template<> inline mpfr::mpreal pi_twice_upper<mpfr::mpreal>() {
			return mpfr::const_pi (mpfr::mpreal::default_prec, MPFR_RNDU)  * 2.0;
		}
	}
	}
  }
}


#endif /* MPREAL_INTERVAL_HPP_ */
