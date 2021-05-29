#ifndef BLAS
#define BLAS

#define BOOST_UBLAS_MOVE_SEMANTICS

#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_sparse.hpp>


#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/triangular.hpp> 
#include <boost/numeric/ublas/lu.hpp>

#include "Util.hpp"

typedef boost::numeric::ublas::compressed_matrix<long double,boost::numeric::ublas::row_major> PrecMatrix;



typedef boost::numeric::ublas::compressed_matrix<double,boost::numeric::ublas::column_major> ColMatrix;



typedef boost::numeric::ublas::mapped_matrix<double,boost::numeric::ublas::row_major> MappedMatrix;

// typedef boost::numeric::ublas::compressed_vector<double> SparseVector;
typedef boost::numeric::ublas::compressed_vector<long double> PrecSparseVector;
typedef boost::numeric::ublas::vector<double> DenseVector;
typedef boost::numeric::ublas::vector<long double> PrecDenseVector;


/************ Memory counting routines ************/
/* count the number of bytes used by a vector     */



/**
 Put the matrix into row-compressed form (eleminate zero rows)
   <M2,nz> := fold(M1)
 */
template<class Matrix>
void fold(Matrix& M1, Matrix& M2, std::vector<unsigned>& M2_nz) {
  typedef typename Matrix::const_iterator1 MatrixConstIterator1;		

  M2_nz.clear();

  // go through the matrix and store the non-zero rows
  for(MatrixConstIterator1 it = M1.begin1(); it!= M1.end1(); ++it) {
    
    if(it.begin() != it.end()) {
      M2_nz.push_back(it.index1());
    }
  }
	
  unsigned nr_of_nz_rows = M2_nz.size();
  
  M2.resize(nr_of_nz_rows,M1.size1(),false);
  
  // build up matrix M2 putting in exactly these non-zero rows
  for(unsigned i=0; i<nr_of_nz_rows; ++i) {
    row(M2,i) = row(M1,M2_nz[i]);
  }
}

template<class Matrix>
void unfold(const Matrix& M1, const std::vector<unsigned>& M1_nz, Matrix& M2) {
	typedef typename Matrix::const_iterator1 MatrixConstIterator1;		

	// unfold M2 to full n x n dimension
	M2.resize(M1.size2(),M1.size2(),false);
	unsigned nr_of_nz_rows = M1_nz.size();
	
    // build up matrix M2 putting in the non-zero rows
	for(unsigned i=0; i<nr_of_nz_rows; ++i) {
		row(M2,M1_nz[i]) = row(M1,i);
    }
}


template<class Vector, class Matrix>
void row_prod(Vector& result, const Vector& v,const Matrix& m) {
  typedef typename Vector::const_iterator VectorConstIterator;		
 for(VectorConstIterator it=v.begin(); it!=v.end(); ++it) {
   const Vector& r = row(m,it.index());
   const long double d(*it);
   for(VectorConstIterator it2=r.begin(); it2!=r.end(); ++it2) {
     result(it2.index()) += d * (*it2);
   }	
 }
}

/*
 * m * v
 */
template<class Vector, class Matrix>
void col_prod(Vector& result, const Vector& v,const Matrix& m, const std::vector<unsigned>& nz) {
  typedef typename Vector::const_iterator VectorConstIterator;
  typedef typename Matrix::const_iterator1 MatrixConstIterator1;		

  result.clear();

  static Vector product(m.size1());
  product.resize(m.size1(),false);
  product.clear();
  assert(m.size2() == v.size());
  product = prec_prod(m,v);
  for(VectorConstIterator it=product.begin(); it!=product.end(); ++it) {
	result(nz[it.index()]) = *it;
  }

}

/*
 * v * m 
 */
template<class Vector, class Matrix>
void row_prod(Vector& result, const Vector& v,const Matrix& m, const std::vector<unsigned>& nz) {
 typedef typename boost::numeric::ublas::matrix_row<Matrix> MatrixRow;
 typedef typename MatrixRow::const_iterator MatrixRowIterator; 	
 typedef typename Vector::const_iterator VectorConstIterator;		

 result.clear(); // first set everythiing to zero
 
 unsigned i = 0;
 for(VectorConstIterator it=v.begin(); it!=v.end(); ++it) {
   unsigned index = it.index();

   // there is no corresponding row, nothing to do
   if(*it == 0 || nz.size () == 0) continue;

   // find the corresponding matrix row
   while( i < nz.size()  && nz[i] < index) {
     ++i;
   } 
   
   // there is no corresponding row, nothing to do
   if( i >= nz.size())
	   break;
   if(nz[i] != index)
	   continue;
   
   MatrixRow r ( row((Matrix&)m,i) );
   const typename Vector::value_type d(*it);
   for(MatrixRowIterator it2=r.begin(); it2!=r.end(); ++it2) {
     result(it2.index()) += d * (*it2);
   }	
 }
}

/*
 * A * B (where B is a row-spare matrix)
 */

template<class Vector, class Matrix>
void row_prod_matrix(Matrix& result, 
		     const Matrix& A,
		     const Matrix& B, 
		     const std::vector<unsigned>& nz) {
 typedef typename boost::numeric::ublas::matrix_row<Matrix> MatrixRow;
 typedef typename MatrixRow::const_iterator MatrixRowIterator; 	
 
 for(unsigned row_index_A=0; row_index_A <A.size(); ++row_index_A) {
   const MatrixRow rowA = row(A,row_index_A);

   unsigned i = 0;
   for(MatrixRowIterator it=rowA.begin(); it!=rowA.end(); ++it) {

     // find the corresponding matrix row
     while( i < nz.size()  && nz[i] < it.index()) {
       ++i;
     } 
   
     // there is no corresponding row, nothing to do
     if(nz[i] != it.index()) continue;
   
     const MatrixRow r = row(B,i);
     const long double d(*it);
     for(MatrixRowIterator it2=r.begin(); it2!=r.end(); ++it2) {
       result(row_index_A,it2.index()) += d * (*it2);
     }	
   }
 }
}

template <class Vector>
void filter(Vector& v, const Precision<typename Vector::value_type>& prec) {
	typedef typename Vector::iterator VectorIterator;
	for(VectorIterator it = v.begin(); it!=v.end();++it)
		*it = fabs(*it) >= prec.absolute_precision ? *it : 0;
}

template<class Matrix>
void filterMatrix(Matrix& M, const Precision<typename Matrix::value_type>& prec) {
	typedef typename Matrix::iterator1 it1_t;
	typedef typename Matrix::iterator2 it2_t;
	for(it1_t it1 = M.begin1(); it1!=M.end1();++it1)
		for(it2_t it2 = it1.begin(); it2!=it1.end();++it2)
			*it2 = fabs(*it2) > prec.absolute_precision ? *it2 : 0;
}

// memory in bytes
template<class Matrix>
long mem(const Matrix& m) {
	return sizeof(int) * m.size1() + m.nnz_capacity() * ( sizeof(double) + sizeof(int)) + sizeof(Matrix);
}

// number of non-zeroes
template<class Matrix>
unsigned nnz(const Matrix& M) {
  typedef typename Matrix::const_iterator1 it1_t;
  typedef typename Matrix::const_iterator2 it2_t;
  unsigned counter = 0;
  for(it1_t it=M.begin1(); it!=M.end1() ; ++it) {
    for(it2_t it2 = it.begin(); it2!=it.end(); ++it2) {
      if(*it2 != 0.0) ++counter;
    }
  }
  return counter;
}

// number of non-zeroes
template<class Matrix>
typename Matrix::value_type col_norm_2(const Matrix& M) {
  typedef typename Matrix::const_iterator1 it1_t;
  typedef typename Matrix::const_iterator2 it2_t;
  typedef typename Matrix::value_type val;
  typedef std::unordered_map<int,val> map_t;


  map_t v;
  for(it1_t it1=M.begin1(); it1!=M.end1() ; ++it1) {
    for(it2_t it2 = it1.begin(); it2!=it1.end(); ++it2) {
      v[it2.index2()] += (*it2) * (*it2);
    }
  }

  val result(0);
  for(typename map_t::const_iterator it=v.begin(); it!=v.end(); ++it) {
	  val entry(sqrt(it->second));
	  result = result < entry ? entry : result;
  }
  return result;
}





double min_nz(const PrecSparseVector& v);

/** Determinant */ 
/*
template<class M> 
double determinant(M const& m) {  
  M mLu(m); 
  boost::numeric::ublas::permutation_matrix<std::size_t> pivots(m.size1()); 
  lu_factorize(mLu, pivots); 
  typename M::value_type det = 1.0; 
  for (std::size_t i=0; i < pivots.size(); ++i) { 
    if (pivots(i) != i) 
      det *= -1.0; 
    det *= mLu(i,i); 
  } 
  return det; 
} 
*/

#endif
