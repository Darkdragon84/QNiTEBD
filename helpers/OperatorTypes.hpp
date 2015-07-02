#ifndef OP_TYPES_H_
#define OP_TYPES_H_

#include <map>
#include <cmath>
#include <iostream>
#include <assert.h>

#include "arma_typedefs.h"

using namespace std;
using std::abs;

template<typename VT>
class
SparseOperator : public SpMat<VT>
{
public:
    SparseOperator(uint d, uint NSites):SpMat<VT>((uint)pow(d,NSites),(uint)pow(d,NSites)),d_(d),NSites_(NSites){}
    SparseOperator(uint d, uint NSites, const SpMat<VT>& mat):SpMat<VT>(mat),d_(d),NSites_(NSites){}
    SparseOperator(uint d, uint NSites, SpMat<VT>&& mat):SpMat<VT>(std::move(mat)),d_(d),NSites_(NSites){}
    SparseOperator(uint d, uint NSites, const umat& locations, const Col<VT>& values):SpMat<VT>(locations,values,(uint)pow(d,NSites),(uint)pow(d,NSites)),d_(d),NSites_(NSites){}

    inline uint GetLocalDim() const {return d_;}
    inline uint GetNSites() const {return NSites_;}
    inline uint GetDim() const {return (uint)pow(d_,NSites_);}

    inline SparseOperator<VT>& operator+=(const SparseOperator<VT>& other);
    inline SparseOperator<VT>& operator-=(const SparseOperator<VT>& other);
    inline SparseOperator<VT>& operator*=(const SparseOperator<VT>& other);
    inline SparseOperator<VT>& operator*=(VT scalar);

    inline SparseOperator<VT> t() const {return SparseOperator<VT>(this->GetLocalDim(),this->GetNSites(),this->SpMat<VT>::t());}

    friend inline SparseOperator<VT> operator*(VT scalar, const SparseOperator<VT>& in) {SparseOperator<VT> out(in);return out*=scalar;}
    friend inline SparseOperator<VT> operator*(const SparseOperator<VT>& in, VT scalar) {SparseOperator<VT> out(in);return out*=scalar;}

    friend inline SparseOperator<VT> operator*(const SparseOperator<VT>& lhs,const SparseOperator<VT>& rhs) {SparseOperator<VT> out(lhs);return out*=rhs;}
    friend inline SparseOperator<VT> operator+(const SparseOperator<VT>& lhs,const SparseOperator<VT>& rhs){SparseOperator<VT> out(lhs);return out+=rhs;}
    friend inline SparseOperator<VT> operator-(const SparseOperator<VT>& lhs,const SparseOperator<VT>& rhs){SparseOperator<VT> out(lhs);return out-=rhs;}
    friend inline SparseOperator<VT> operator-(const SparseOperator<VT>& in){SparseOperator<VT> out(in);return out*=(-1.);}

    static inline SparseOperator<VT> SpId(uint d, uint N) {return SparseOperator<VT>(d,N,speye((uint)pow(d,N),(uint)pow(d,N)));}

protected:
    uint d_,NSites_;
};

typedef SparseOperator<Real> RSpOp;
typedef SparseOperator<Complex> CSpOp;

template<typename VT>
inline
SparseOperator<VT>&
SparseOperator<VT>::operator+=(const SparseOperator<VT>& other)
{
    assert(this->GetLocalDim() == other.GetLocalDim() && this->GetNSites() == other.GetNSites());
    this->SpMat<VT>::operator+=(other);
    return *this;
}

template<typename VT>
inline
SparseOperator<VT>&
SparseOperator<VT>::operator-=(const SparseOperator<VT>& other)
{
    assert(this->GetLocalDim() == other.GetLocalDim() && this->GetNSites() == other.GetNSites());
    this->SpMat<VT>::operator-=(other);
    return *this;
}

template<typename VT>
inline
SparseOperator<VT>&
SparseOperator<VT>::operator*=(const SparseOperator<VT>& other)
{
    assert(this->GetLocalDim() == other.GetLocalDim() && this->GetNSites() == other.GetNSites());
    this->SpMat<VT>::operator*=(other);
    return *this;
}

template<typename VT>
inline
SparseOperator<VT>&
SparseOperator<VT>::operator*=(VT scalar)
{
    this->SpMat<VT>::operator*=(scalar);
    return *this;
}

/**< external functuins ****************************************************************** */

template<typename VT>
SparseOperator<VT>
kron(const SparseOperator<VT>& lhs, const SparseOperator<VT>& rhs)
{
    assert(lhs.GetLocalDim() == rhs.GetLocalDim());
    uint rdim = rhs.GetDim();
    uint newdim = lhs.n_nonzero * rhs.n_nonzero;
    umat locs(2,newdim);
    Col<VT> vals(newdim);

    uint ct=0;
    for (typename SparseOperator<VT>::const_iterator lhit=lhs.begin(); lhit!=lhs.end(); ++lhit)
    {
        for (typename SparseOperator<VT>::const_iterator rhit=rhs.begin(); rhit!=rhs.end(); ++rhit)
        {
            locs(0,ct)=lhit.row()*rdim + rhit.row();
            locs(1,ct)=lhit.col()*rdim + rhit.col();
            vals(ct)=(*lhit)*(*rhit);
            ++ct;
        }
    }
    return SparseOperator<VT>(lhs.GetLocalDim(),lhs.GetNSites() + rhs.GetNSites(),locs,vals);
}

template<typename VT>
inline
SparseOperator<VT>
expmat(const SparseOperator<VT>& in)
{
    return SparseOperator<VT>(in.GetLocalDim(),in.GetNSites(),SpMat<VT>(expmat(Mat<VT>(in))));
}

/**< SOME SPECIAL OPERATORS *************************************************************** */



#endif //OP_TYPES_H_
