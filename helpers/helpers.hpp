#ifndef HELPERS_H_
#define HELPERS_H_

#include <cmath>
#include <iostream>
#include <assert.h>
//#include <utility>
//#include <time.h>

#include "arma_typedefs.h"

//#define pseudoinv(A,B) (abs(A) < (B)  ? (0) : (1./A))

using std::cout;
using std::endl;


/**< DENSE MATRIX MULTIPLICATIONS AND DIVISIONS BY SCHMIDT VALUES *********************************************************/
/**< multiplication */
template<typename T>
inline
void
MultMatLamRight(Mat<T>& in, const RVecType& lam)
{
    assert(in.n_cols==lam.n_elem);
//    for (uint i=0;i<lam.n_elem;++i) in.unsafe_col(i)*=lam(i);
    for (uint i=0;i<lam.n_elem;++i) in.col(i)*=lam(i);
}

template<typename T>
inline
void
MultMatLamLeft(const RVecType& lam, Mat<T>& in)
{
    assert(in.n_rows==lam.n_elem);
    for (uint i=0;i<lam.n_elem;++i) in.row(i)*=lam(i);
}

/**< division */
template<typename T>
inline
void
DivMatLamRight(Mat<T>& in, const RVecType& lam, double thresh=1e-14)
{
    assert(in.n_cols==lam.n_elem);
    for (uint i=0;i<lam.n_elem;++i)
    {
        if(lam(i)<thresh)cerr<<"DivMatLamRight: division by "<<lam(i)<<endl;
//        in.unsafe_col(i)/=lam(i);
        in.col(i)/=lam(i);
    }
}

template<typename T>
inline
void
DivMatLamLeft(const RVecType& lam, Mat<T>& in, double thresh=1e-14)
{
    assert(in.n_rows==lam.n_elem);
    for (uint i=0;i<lam.n_elem;++i)
    {
        if(lam(i)<thresh)cerr<<"DivMatLamLeft: division by "<<lam(i)<<endl;
        in.row(i)/=lam(i);
    }
}

/**< modifying operations */
/**< multiplication */
template<typename T>
inline
Mat<T>&
operator<(Mat<T>& in, const RVecType& lam) {MultMatLamRight(in,lam); return in;}

template<typename T>
inline
Mat<T>&
operator>(const RVecType& lam, Mat<T>& in) {MultMatLamLeft(lam,in); return in;}

/**< division */
template<typename T>
inline
Mat<T>&
operator>(Mat<T>& in, const RVecType& lam) {DivMatLamRight(in,lam); return in;}

template<typename T>
inline
Mat<T>&
operator<(const RVecType& lam, Mat<T>& in) {DivMatLamLeft(lam,in); return in;}

/**< non-modifying operations (create copies first)*/
/**< multiplication */
template<typename T>
inline
Mat<T>
operator<<(const Mat<T>& in, const RVecType& lam) {Mat<T> out(in); MultMatLamRight(out,lam); return out;}

template<typename T>
inline
Mat<T>
operator>>(const RVecType& lam, const Mat<T>& in) {Mat<T> out(in); MultMatLamLeft(lam,out); return out;}

/**< division */
template<typename T>
inline
Mat<T>
operator>>(const Mat<T>& in, const RVecType& lam) {Mat<T> out(in); DivMatLamRight(out,lam); return out;}

template<typename T>
inline
Mat<T>
operator<<(const RVecType& lam, const Mat<T>& in) {Mat<T> out(in); DivMatLamLeft(lam,out); return out;}

/**< dot product for matrices */
//template<typename T>
//inline
//T
//dot(const Mat<T>& lhs, const Mat<T>& rhs)
//{
//    assert(lhs.n_elem==rhs.n_elem);
//    return dot(Col<T>(lhs.memptr(),lhs.n_elem),Col<T>(rhs.memptr(),rhs.n_elem));
//}
//
//template<typename T>
//inline
//Real
//dot(const Mat<T>& mat){return real(dot(mat,mat));}

#endif
