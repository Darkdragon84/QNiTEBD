#ifndef BLOCK_OBJ_H
#define BLOCK_OBJ_H

#include <map>
#include "helpers.hpp"

using std::ostream;
using std::endl;
using std::cout;

/**< BLOCK OBJECTS AS BUILDING BLOCKS FOR MPS  */

/**< Block Schmidt Value Type *****************************************************************************************
defined as a map of Real valued vectors (RVecType), where the quantum number(s) of the blocks are the keys for the map
*/
template<typename KT>
class BlockLambda : public std::map<KT,RVecType>
{
public:
    void ShowDims(const std::string& name="") const;
    inline uint GetTotalSize() const;
};

template<typename KT>
void
BlockLambda<KT>::ShowDims(const std::string& name) const
{
    if (name!="")cout<<name<<endl;
    if (!this->empty())
    {
        for (const auto& it : *this) cout<<it.first<<": "<<it.second.size()<<endl;
    }
    else cout<<"---"<<endl;
}

template<typename KT>
inline
uint
BlockLambda<KT>::GetTotalSize() const
{
    uint dim=0;
    for (const auto it : *this) dim+=it.second.n_elem;
    return dim;
}

template<typename KT>
ostream& operator<<(ostream& os, const BlockLambda<KT>& lam)
{
    if (!lam.empty())
    {
        for (const auto& it : lam)
        {
            os<<it.first<<endl;
            os<<it.second;
        }
    }
    else os<<"---"<<endl;
    return os;
}

template<typename KT>
inline
BlockLambda<KT>&
operator*=(BlockLambda<KT>& lam, Real val)
{
    for (auto& it : lam) it.second*=val;
    return lam;
}

template<typename KT>
inline
BlockLambda<KT>
operator*(Real val, const BlockLambda<KT>& lamin)
{
    BlockLambda<KT> lamout(lamin);
    return lamout*=val;
}

template<typename KT>
inline
BlockLambda<KT>
operator*(const BlockLambda<KT>& lamin, Real val)
{
    BlockLambda<KT> lamout(lamin);
    return lamout*=val;
}

template<typename KT>
inline
BlockLambda<KT>&
operator/=(BlockLambda<KT>& lam, Real val)
{
    for (auto& it : lam) it.second/=val;
    return lam;
}

template<typename KT>
inline
BlockLambda<KT>
operator/(const BlockLambda<KT>& lamin, Real val)
{
    BlockLambda<KT> lamout(lamin);
    return lamout/=val;
}


template<typename KT>
inline
Real
norm(const BlockLambda<KT>& lam)
{
    Real nrm2=0;
    for (const auto& it : lam) nrm2+=dot(it.second,it.second);
    return sqrt(nrm2);
}

template<typename KT>
inline
BlockLambda<KT>
sqrt(const BlockLambda<KT>& in)
{
    BlockLambda<KT> out;
    for (const auto& it : in) out.emplace_hint(out.end(),it.first,sqrt(it.second));
    return out;
}

template<typename KT,typename VT>
inline
BlockLambda<KT>
pow(const BlockLambda<KT>& in, VT expo)
{
    BlockLambda<KT> out;
    for (const auto& it : in) out.emplace_hint(out.end(),it.first,pow(it.second,expo));
    return out;
}

/** BlockDiagonal Matrix Type ************************************************************************
 *  the ingoing and outgoing quantum number is the same, therefore the second needs not to be stored
 */

template<typename KT,typename VT>
class BlockDiagMat : public std::map<KT,Mat<VT> >
{
    typedef typename std::map<KT,Mat<VT> > maptype;
public:
    BlockDiagMat() = default;
    BlockDiagMat(const BlockDiagMat& in) = default;
    BlockDiagMat(BlockDiagMat&& in) = default;
    BlockDiagMat& operator=(const BlockDiagMat& other) & = default;
    BlockDiagMat& operator=(BlockDiagMat&& other) & = default;
    BlockDiagMat(const BlockLambda<KT>& lam);
    void ShowDims(const std::string& name="") const;

    inline BlockDiagMat& operator*=(VT scalar);
    inline BlockDiagMat& operator/=(VT scalar);

    friend BlockDiagMat operator*(VT scalar, const BlockDiagMat& in) {BlockDiagMat out(in);return out*=scalar;}
    friend BlockDiagMat operator*(const BlockDiagMat& in, VT scalar) {BlockDiagMat out(in);return out*=scalar;}
};

template<typename KT,typename VT>
BlockDiagMat<KT,VT>::BlockDiagMat(const BlockLambda<KT>& lam)
{
    for (const auto& it : lam) this->emplace_hint(this->end(),it.first,diagmat(it.second));
}

/**< helper function to show all present quantum number sectors */
template<typename KT,typename VT>
void
BlockDiagMat<KT,VT>::ShowDims(const std::string& name) const
{
    if (name!="")cout<<name<<endl;
    if (!this->empty())
    {
        for (const auto& it : *this) cout<<it.first<<": "<<it.second.n_rows<<"x"<<it.second.n_cols<<endl;
    }
    else cout<<"---"<<endl;
}

template<typename KT,typename VT>
BlockDiagMat<KT,VT>&
BlockDiagMat<KT,VT>::operator*=(VT scalar)
{
    for (auto& it : *this) it.second*=scalar;
    return *this;
}

template<typename KT,typename VT>
BlockDiagMat<KT,VT>&
BlockDiagMat<KT,VT>::operator/=(VT scalar)
{
    for (auto& it : *this) it.second/=scalar;
    return *this;
}

/**< output stream operator overload for screen output */
template<typename KT,typename VT>
ostream&
operator<<(ostream& os, const BlockDiagMat<KT,VT>& M)
{
    if (!M.empty())
    {
        for (const auto& it : M)
        {
            os<<it.first<<":"<<endl;
            os<<it.second<<endl;
        }
    }
    else os<<"---"<<endl;
    return os;
}

/**< TRIVIAL ARITHMETICS FOR BLOCKDIAGMATS ***********************************************************/

template<typename KT,typename VT>
inline
BlockDiagMat<KT,VT>
trans(const BlockDiagMat<KT,VT>& in)
{
    BlockDiagMat<KT,VT> out(in);
    for (auto& it : out) it.second.t();
    return out;
}

template<typename KT,typename VT>
inline
Real
norm(const BlockDiagMat<KT,VT>& in)
{
    Real nrm2=0;
    for (const auto& it : in) nrm2+=abs(dot(it.second,it.second));
//    for (const auto& it : in) nrm2+=pow(norm(it.second,"fro"),2);
    return sqrt(nrm2);
}

template<typename KT,typename VT>
inline
BlockDiagMat<KT,VT>
operator*(VT scalar, const BlockDiagMat<KT,VT>& in)
{
    BlockDiagMat<KT,VT> out;
    if (!in.empty())
    {
        for (const auto& init : in) out.emplace_hint(out.end(),init.first,scalar*init.second);
    }
    return out;
}

template<typename KT,typename VT>
inline
BlockDiagMat<KT,VT>
operator+(const BlockDiagMat<KT,VT>& lhs, const BlockDiagMat<KT,VT>& rhs)
{
    BlockDiagMat<KT,VT> out;
    for (const auto& lhit : lhs)
    {
        const auto rhit = rhs.find(lhit.first);
        if (rhit!=rhs.end()) out.emplace_hint(out.end(),lhit.first,lhit.second + rhit->second);
    }
    return out;
}

template<typename KT,typename VT>
inline
BlockDiagMat<KT,VT>
operator-(const BlockDiagMat<KT,VT>& lhs, const BlockDiagMat<KT,VT>& rhs)
{
    BlockDiagMat<KT,VT> out;
    for (const auto& lhit : lhs)
    {
        const auto rhit = rhs.find(lhit.first);
        if (rhit!=rhs.end()) out.emplace_hint(out.end(),lhit.first,lhit.second - rhit->second);
    }
    return out;
}

template<typename KT,typename VT>
inline
VT
trace(const BlockDiagMat<KT,VT>& mat)
{
    VT out=0;
    for (const auto& it : mat) out+=trace(it.second);
    return out;
}

/**< BASIC MULTIPLICATION AND DIVISION BY DIAGONAL MATRICES (E.G. SCHMIDT VALUES) ******************************************************/

template<typename KT,typename VT>
BlockDiagMat<KT,VT>&
MultBlockDiagMatLamLeft(const BlockLambda<KT>& lam, BlockDiagMat<KT,VT>& in)
{
    if (!lam.empty())
    {
        for (auto& matit : in)
        {
            const auto lamit = lam.find(matit.first);
            if (lamit != lam.end()) MultMatLamLeft(lamit->second,matit.second);
        }
    }
    return in;
}


template<typename KT,typename VT>
BlockDiagMat<KT,VT>&
MultBlockDiagMatLamRight(BlockDiagMat<KT,VT>& in, const BlockLambda<KT>& lam)
{
    if (!lam.empty())
    {
        for (auto& matit : in)
        {
            const auto lamit = lam.find(matit.first);
            if (lamit!=lam.end()) matit.second < (lamit->second);
        }
    }
    return in;
}

template<typename KT,typename VT>
BlockDiagMat<KT,VT>&
DivBlockDiagMatLamLeft(const BlockLambda<KT>& lam,BlockDiagMat<KT,VT>& in)
{
    if (!lam.empty())
    {
        for (auto& matit : in)
        {
            const auto lamit = lam.find(matit.first);
            if (lamit!=lam.end()) (lamit->second) < matit.second;
        }
    }
    return in;
}


template<typename KT,typename VT>
BlockDiagMat<KT,VT>&
DivBlockDiagMatLamRight(BlockDiagMat<KT,VT>& in, const BlockLambda<KT>& lam)
{
    if (!lam.empty())
    {
        for (auto& matit : in)
        {
            const auto lamit = lam.find(matit.first);
            if (lamit!=lam.end()) matit.second > (lamit->second);
        }
    }
    return in;
}

/// modifying
template<typename KT,typename VT>
inline BlockDiagMat<KT,VT>& operator>(const BlockLambda<KT>& lam, BlockDiagMat<KT,VT>& in) {return MultBlockDiagMatLamLeft(lam,in);}

template<typename KT,typename VT>
inline BlockDiagMat<KT,VT>& operator<(BlockDiagMat<KT,VT>& in, const BlockLambda<KT>& lam) {return MultBlockDiagMatLamRight(in,lam);}

template<typename KT,typename VT>
inline BlockDiagMat<KT,VT>& operator<(const BlockLambda<KT>& lam, BlockDiagMat<KT,VT>& in) {return DivBlockDiagMatLamLeft(lam,in);}

template<typename KT,typename VT>
inline BlockDiagMat<KT,VT>& operator>(BlockDiagMat<KT,VT>& in, const BlockLambda<KT>& lam) {return DivBlockDiagMatLamRight(in,lam);}

/// non-modifying
template<typename KT,typename VT>
inline BlockDiagMat<KT,VT> operator>>(const BlockLambda<KT>& lam, const BlockDiagMat<KT,VT>& in) {BlockDiagMat<KT,VT> out(in); return MultBlockDiagMatLamLeft(lam,out);}

template<typename KT,typename VT>
inline BlockDiagMat<KT,VT> operator<<(const BlockDiagMat<KT,VT>& in, const BlockLambda<KT>& lam) {BlockDiagMat<KT,VT> out(in); return MultBlockDiagMatLamRight(out,lam);}

template<typename KT,typename VT>
inline BlockDiagMat<KT,VT> operator<<(const BlockLambda<KT>& lam, const BlockDiagMat<KT,VT>& in) {BlockDiagMat<KT,VT> out(in); return DivBlockDiagMatLamLeft(lam,out);}

template<typename KT,typename VT>
inline BlockDiagMat<KT,VT> operator>>(const BlockDiagMat<KT,VT>& in, const BlockLambda<KT>& lam) {BlockDiagMat<KT,VT> out(in); return DivBlockDiagMatLamRight(out,lam);}


/** Block Matrix Type ************************************************************************
 *  Defined as a map of matrices of numerical type VT. Each block has an ingoing and outgoing (set of)
 *  quantum number(s) of type KT. The blocks are ordered according to the ingoing quantum number (serves
 *  as map key), the mapped type is then a pair consisting of the outgoing quantum number(s) and the matrix
 */
template<typename KT,typename VT>
using QMatPair = std::pair<KT,Mat<VT> >;

template<typename KT,typename VT>
class BlockMat : public std::map<KT,QMatPair<KT,VT> >
{
public:
    typedef typename std::map<KT,QMatPair<KT,VT> > maptype;
    typedef typename std::map<KT,QMatPair<KT,VT> >::mapped_type entry;
    typedef typename std::map<KT,QMatPair<KT,VT> >::value_type value;

//    BlockMat(BlockMat&& in):maptype(std::move(in)) {DOUT("BlockMat copied");};
//    BlockMat(const BlockMat& in):maptype(in) {DOUT("BlockMat copied");};
    void ShowDims(const std::string& name="") const;
};


/**< ACCESS HELPERS FOR BlockMat -------------------- */
/**< if O is a Block in a BlockMat object, call like Qin(O), Qout(O), QMat(O) to get the respective quantum numbers */
template<typename KT,typename VT>
inline
const KT&
Qin(const std::pair<const KT,QMatPair<KT,VT> >& E)
{
    return E.first;
}

template<typename KT,typename VT>
inline
const KT&
Qout(const std::pair<const KT,QMatPair<KT,VT> >& E)
{
    return E.second.first;
}

template<typename KT,typename VT>
inline
Mat<VT>&
QMat(std::pair<const KT,std::pair<KT,Mat<VT> > >& E)
{
    return E.second.second;
}

template<typename KT,typename VT>
inline
const Mat<VT>&
QMat(const std::pair<const KT,std::pair<KT,Mat<VT> > >& E)
{
    return E.second.second;
}

/**< helper function to show all present quantum number sectors */
template<typename KT,typename VT>
void
BlockMat<KT,VT>::ShowDims(const std::string& name) const
{
    if (name!="")cout<<name<<endl;
    if (!this->empty())
    {
        for (const auto& it : *this)cout<<"("<<Qin(it)<<","<<Qout(it)<<"): "<<QMat(it).n_rows<<"x"<<QMat(it).n_cols<<endl;
    }
    else cout<<"---"<<endl;
}

/**< output stream operator overload for screen output */
template<typename KT,typename VT>
ostream&
operator<<(ostream& os, const BlockMat<KT,VT>& M)
{
    if (!M.empty())
    {
        for (const auto& it : M)
        {
            os<<"("<<Qin(it)<<","<<Qout(it)<<"):"<<endl;
            os<<QMat(it)<<endl;
        }
    }
    else os<<"---"<<endl;
    return os;
}

/**< ADDITION AND SUBTRACTION FOR BLOCK MATRICES */
//template<typename KT,typename VT>
//BlockMat<KT,VT>&
//operator+=(BlockMat<KT,VT>& lhs, const BlockMat<KT,VT>& rhs)
//{
//    for (const auto& rit : rhs)
//    {
//        auto lit = lhs.find(Qin(rit));
//        if (lit==lhs.end()) lhs.emplace(Qin(rit),QMatPair(Qout(rit),QMat(rit)));
//        else (*lit)+=
//    }
//}

/**< MULTIPLICATION ROUTINES FOR BLOCK MATRICES ****************************************************/

/**< multiplication operator for two BlockMats */
template<typename KT,typename VT>
BlockMat<KT,VT>
operator*(const BlockMat<KT,VT>& lhs, const BlockMat<KT,VT>& rhs)
{
    BlockMat<KT,VT> out;
    for (const auto& lit : lhs)
    {
        const auto rit = rhs.find(Qout(lit));
        if (rit!=rhs.end()) out.emplace_hint(out.end(),Qin(lit),QMatPair<KT,VT>(Qout(*rit) , QMat(lit) * QMat(*rit) ));
    }
    return out;
}

/**< multiplication operator for BlockMat (left) and BlockDiagMat (right) */
template<typename KT,typename VT1,typename VT2> /// it could in principle happen that BlockMat is complex and BlockDiagMat is real
BlockMat<KT,VT1>
operator*(const BlockMat<KT,VT1>& lhs, const BlockDiagMat<KT,VT2>& rhs)
{
    BlockMat<KT,VT1> out;
    for (const auto& lit : lhs)
    {
        const auto rit = rhs.find(Qout(lit));
        if (rit!=rhs.end()) out.emplace_hint(out.end(),Qin(lit),QMatPair<KT,VT1>(Qout(lit) , QMat(lit) * rit->second ));
    }
    return out;
}

/**< multiplication operator for BlockDiagMat (left) and BlockMat (right) */
template<typename KT,typename VT1,typename VT2> /// it could in principle happen that BlockMat is complex and BlockDiagMat is real
BlockMat<KT,VT1>
operator*(const BlockDiagMat<KT,VT2>& lhs, const BlockMat<KT,VT1>& rhs)
{
    BlockMat<KT,VT1> out;
    for (const auto& rit : rhs)
    {
        const auto lit = lhs.find(Qin(rit));
        if (lit!=lhs.end()) out.emplace_hint(out.end(),Qin(rit),QMatPair<KT,VT1>(Qout(rit) , lit->second * QMat(rit) ));
    }
    return out;
}

/**< multiplication operator for two BlockDiagMats */
template<typename KT,typename VT>
BlockDiagMat<KT,VT>
operator*(const BlockDiagMat<KT,VT>& lhs, const BlockDiagMat<KT,VT>& rhs)
{
    BlockDiagMat<KT,VT> out;
    for (const auto& lit : lhs)
    {
        const auto rit = rhs.find(lit.first);
        if (rit!=rhs.end()) out.emplace_hint(out.end(),lit.first, lit.second*rit->second );
    }
    return out;
}



/**< BASIC MULTIPLICATION AND DIVISION BY DIAGONAL MATRICES (E.G. SCHMIDT VALUES) ******************************************************/

template<typename KT,typename VT>
BlockMat<KT,VT>&
MultBlockMatLamLeft(const BlockLambda<KT>& lam,BlockMat<KT,VT>& in)
{
    if (!lam.empty())
    {
        for (auto& matit : in)
        {
            const auto lamit = lam.find(Qin(matit));
            if (lamit!=lam.end()) (lamit->second) > QMat(matit);
        }
    }
    return in;
}


template<typename KT,typename VT>
BlockMat<KT,VT>&
MultBlockMatLamRight(BlockMat<KT,VT>& in, const BlockLambda<KT>& lam)
{
    if (!lam.empty())
    {
        for (auto& matit : in)
        {
            const auto lamit = lam.find(Qout(matit));
            if (lamit!=lam.end()) QMat(matit) < (lamit->second);
        }
    }
    return in;
}

template<typename KT,typename VT>
BlockMat<KT,VT>&
DivBlockMatLamLeft(const BlockLambda<KT>& lam,BlockMat<KT,VT>& in)
{
    if (!lam.empty())
    {
        for (auto& matit : in)
        {
            const auto lamit = lam.find(Qin(matit));
            if (lamit!=lam.end()) (lamit->second) < QMat(matit);
        }
    }
    return in;
}


template<typename KT,typename VT>
BlockMat<KT,VT>&
DivBlockMatLamRight(BlockMat<KT,VT>& in, const BlockLambda<KT>& lam)
{
    if (!lam.empty())
    {
        for (auto& matit : in)
        {
            const auto lamit = lam.find(Qout(matit));
            if (lamit!=lam.end()) QMat(matit) > (lamit->second);
        }
    }
    return in;
}

/// modifying
template<typename KT,typename VT>
inline BlockMat<KT,VT>& operator>(const BlockLambda<KT>& lam, BlockMat<KT,VT>& in) {return MultBlockMatLamLeft(lam,in);}

template<typename KT,typename VT>
inline BlockMat<KT,VT>& operator<(BlockMat<KT,VT>& in, const BlockLambda<KT>& lam) {return MultBlockMatLamRight(in,lam);}

template<typename KT,typename VT>
inline BlockMat<KT,VT>& operator<(const BlockLambda<KT>& lam, BlockMat<KT,VT>& in) {return DivBlockMatLamLeft(lam,in);}

template<typename KT,typename VT>
inline BlockMat<KT,VT>& operator>(BlockMat<KT,VT>& in, const BlockLambda<KT>& lam) {return DivBlockMatLamRight(in,lam);}

/// non-modifying
template<typename KT,typename VT>
inline BlockMat<KT,VT> operator>>(const BlockLambda<KT>& lam, const BlockMat<KT,VT>& in) {BlockMat<KT,VT> out(in); return MultBlockMatLamLeft(lam,out);}

template<typename KT,typename VT>
inline BlockMat<KT,VT> operator<<(const BlockMat<KT,VT>& in, const BlockLambda<KT>& lam) {BlockMat<KT,VT> out(in); return MultBlockMatLamRight(out,lam);}

template<typename KT,typename VT>
inline BlockMat<KT,VT> operator<<(const BlockLambda<KT>& lam, const BlockMat<KT,VT>& in) {BlockMat<KT,VT> out(in); return DivBlockMatLamLeft(lam,out);}

template<typename KT,typename VT>
inline BlockMat<KT,VT> operator>>(const BlockMat<KT,VT>& in, const BlockLambda<KT>& lam) {BlockMat<KT,VT> out(in); return DivBlockMatLamRight(out,lam);}

#endif // BLOCK_OBJ_H
