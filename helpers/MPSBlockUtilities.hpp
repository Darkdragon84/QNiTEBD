//#ifndef _USE_SYMMETRIES_
//#define _USE_SYMMETRIES_
//#endif // _USE_SYMMETRIES_

#ifndef MPS_BLOCK_UTIL_H
#define MPS_BLOCK_UTIL_H

#include <utility>

#include "MPSBlockMat.hpp"
#include "OperatorTypes.hpp"
#include "helpers.hpp"
#include "IKey.hpp"

using std::cout;
using std::endl;
using std::map;

/**< Apply Operator onto MPS */
template<uint N,typename KT,typename VT>
MPSBlockMat<KT,VT>
ApplyOperator(const MPSBlockMat<KT,VT>& MPSin, const RSpOp& op, const ItoKey<N>& I2K)
{
    assert(MPSin.GetNSites() == N && op.GetNSites()==N);
    assert(MPSin.GetLocalDim() == op.GetLocalDim());
    assert(MPSin.GetLocalDim() == I2K.GetLocalDim());
//    assert(MPSin.GetLocalDim() == op.GetLocalDim() && MPSin.GetNSites() == op.GetNSites());
//    assert(I2K.GetLocalDim() == MPSin.GetLocalDim() && I2K.GetNSites() == MPSin.GetNSites());
    MPSBlockMat<KT,VT> MPSout(MPSin.GetLocalDim(),MPSin.GetNSites());
    uint ii,jj;
    for (RSpOp::const_iterator opit=op.begin(); opit!=op.end(); opit++)
    {
        ii=opit.row();
        jj=opit.col();
        if (I2K[ii] == I2K[jj])
        {
            for (const auto init : MPSin[jj])
            {
//                auto outit = MPSout[ii].find(Qin(init));
//                if (outit==MPSout[ii].end()) MPSout[ii].emplace(Qin(init),QMatPair<KT,VT>(Qout(init), (*opit) * QMat(init)));
//                else QMat(*outit) += (*opit) * QMat(init);
                auto outit = MPSout[ii].lower_bound(Qin(init));
                if (outit==MPSout[ii].end() || outit->first!=Qin(init)) MPSout[ii].emplace_hint(outit,Qin(init),QMatPair<KT,VT>(Qout(init), (*opit) * QMat(init)));
                else QMat(*outit) += (*opit) * QMat(init);
            }
        }
    }
    return MPSout;
}

/**< Regular Transfer Matrix Operations ****************************************************************************************/

/**< Apply regular transfer matrix, generated by A, to the left onto unity (i.e. only \sum_s A_s' * A_s
 *   since there is no input matrix, the  output matrix is a BlockDiagMat */
template<typename KT,typename VT>
BlockDiagMat<KT,VT>
ApplyTMLeft(const MPSBlockMat<KT,VT>& A)
{
    BlockDiagMat<KT,VT> out;
    for (const auto& Ait : A) /// loop through phys. indices of A, Ait has type map of matrices, i.e. BlockMat !!
    {
        if (!Ait.empty()) /// see if there are any entries (i.e. symmetry sectors for current phys. index)
        {
            for (const auto& mapit : Ait) /// loop through symmetry sectors of current phys. index (matit is of type QMatPair)
            {
//                auto outit = out.find(Qout(mapit)); /// look for current symmetry sector (result of A_s'*A_s has OUTgoing symmetry label of A_s)
//                if (outit==out.end()) out.emplace(Qout(mapit), QMat(mapit).t() * QMat(mapit) ); /// if symmetry sector not yet present, create
//                else outit->second += QMat(mapit).t() * QMat(mapit); /// otherwise add result to already present contribution
                auto outit = out.lower_bound(Qout(mapit)); /// look for current symmetry sector (result of A_s'*A_s has OUTgoing symmetry label of A_s)
                if (outit==out.end() || outit->first!=Qout(mapit)) out.emplace_hint(outit,Qout(mapit), QMat(mapit).t() * QMat(mapit) ); /// if symmetry sector not yet present, create
                else outit->second += QMat(mapit).t() * QMat(mapit); /// otherwise add result to already present contribution
            }
        }
    }
    return out;
}

/**< Apply regular transfer matrix, generated by A, to the right onto unity (i.e. only \sum_s A_s * A_s'
 *  since there is no input matrix, the  output matrix is a BlockDiagMat */
template<typename KT,typename VT>
BlockDiagMat<KT,VT>
ApplyTMRight(const MPSBlockMat<KT,VT>& A)
{
    BlockDiagMat<KT,VT> out;
    for (const auto& Ait : A) /// loop through phys. indices of A, Ait has type map of matrices, i.e. BlockMat !!
    {
        if (!Ait.empty()) /// see if there are any entries (i.e. symmetry sectors for current phys. index)
        {
            for (const auto& mapit : Ait) /// loop through symmetry sectors of current phys. index (matit is of type QMatPair)
            {
//                auto outit = out.find(Qin(mapit)); /// look for current symmetry sector (result of A_s * A_s' has INgoing symmetry label of A_s)
//                if (outit==out.end()) out.emplace(Qin(mapit), QMat(mapit) * QMat(mapit).t() ); /// if symmetry sector not yet present, create
//                else outit->second += QMat(mapit) * QMat(mapit).t(); /// otherwise add result to already present contribution
                auto outit = out.lower_bound(Qin(mapit)); /// look for current symmetry sector (result of A_s * A_s' has INgoing symmetry label of A_s)
                if (outit==out.end() || outit->first!=Qin(mapit)) out.emplace_hint(outit,Qin(mapit), QMat(mapit) * QMat(mapit).t() ); /// if symmetry sector not yet present, create
                else outit->second += QMat(mapit) * QMat(mapit).t(); /// otherwise add result to already present contribution
            }
        }
    }
    return out;
}

//template<typename KT,typename VT>
//inline
//BlockDiagMat<KT,VT>
//ApplyTMLeft(const MPSBlockMat<KT,VT>& A, const BlockDiagMat<KT,VT>& in){return ApplyTMmixedLeft(in*A,A);}
template<typename KT,typename VTMPS,typename VTmat>
BlockDiagMat<KT,VTmat>
ApplyTMLeft(const MPSBlockMat<KT,VTMPS>& A, const BlockDiagMat<KT,VTmat>& in)
{
    BlockDiagMat<KT,VTmat> out;
    for (const auto& Ait : A)
    {
        if (!Ait.empty()) /// see if there are any symmetry sectors for current physical index s
        {
            for (const auto& matit : Ait) /// loop through symmetry sectors of A_s
            {
                const auto init = in.find(Qin(matit)); /// see if current INgoing symmetry sector is present in input matrix
                if (init!=in.end()) /// only do something if found
                {
//                    auto outit = out.find(Qout(matit)); /// see if current OUTgoing symmetry sector already is present in output matrix
//                    if (outit==out.end()) out.emplace(Qout(matit),QMat(matit).t() * ((init->second) * QMat(matit)) ); /// if not, create and store current result
//                    else outit->second += QMat(matit).t() * ((init->second) * QMat(matit)); /// if yes, add current result to already present contribution
                    auto outit = out.lower_bound(Qout(matit)); /// see if current OUTgoing symmetry sector already is present in output matrix
                    if (outit==out.end() || outit->first!=Qout(matit)) out.emplace_hint(outit,Qout(matit),QMat(matit).t() * ((init->second) * QMat(matit)) ); /// if not, create and store current result
                    else outit->second += QMat(matit).t() * ((init->second) * QMat(matit)); /// if yes, add current result to already present contribution
                }
                else DOUT(Qin(matit)<<" not found"<<endl);
            }
        }
    }
    return out;
}

template<typename KT,typename VTMPS,typename VTmat>
BlockDiagMat<KT,VTmat>
ApplyTMRight(const MPSBlockMat<KT,VTMPS>& A, const BlockDiagMat<KT,VTmat>& in)
{
    BlockDiagMat<KT,VTmat> out;
    for (const auto Ait : A)
    {
        if (!Ait.empty()) /// see if there are any symmetry sectors for current physical index s
        {
            for (const auto& matit : Ait) /// loop through symmetry sectors of A_s
            {
                const auto init = in.find(Qout(matit)); /// see if current OUTgoing symmetry sector is present in input matrix
                if (init!=in.end()) /// only do something if found
                {
//                    auto outit = out.find(Qin(matit)); /// see if current INgoing symmetry sector already is present in output matrix
//                    if (outit==out.end()) out.emplace(Qin(matit), (QMat(matit) * init->second) * QMat(matit).t() ); /// if not, create and store current result
//                    else outit->second += (QMat(matit) * init->second) * QMat(matit).t(); /// if yes, add current result to already present contribution
                    auto outit = out.lower_bound(Qin(matit)); /// see if current INgoing symmetry sector already is present in output matrix
                    if (outit==out.end() || outit->first!=Qin(matit)) out.emplace_hint(outit,Qin(matit), (QMat(matit) * init->second) * QMat(matit).t() ); /// if not, create and store current result
                    else outit->second += (QMat(matit) * init->second) * QMat(matit).t(); /// if yes, add current result to already present contribution
                }
            }
        }
    }
    return out;
}

/**< This version is meant to be used in TMEigs only, where the memory of out is preallocated for speed */
template<typename KT,typename VTMPS,typename VTmat>
void
ApplyTMLeft(const MPSBlockMat<KT,VTMPS>& A, const BlockDiagMat<KT,VTmat>& in, BlockDiagMat<KT,VTmat>& out)
{
    /**< ATTENTION: IF OUT IS NOT COMPLETELY EMPTY, MAKE SURE IT IS INITIALIZED TO ZERO */
    for (const auto& Ait : A)
    {
        if (!Ait.empty()) /// see if there are any symmetry sectors for current physical index s
        {
            for (const auto& matit : Ait) /// loop through symmetry sectors of A_s
            {
                const auto init = in.find(Qin(matit)); /// see if current INgoing symmetry sector is present in input matrix
                if (init!=in.end()) /// only do something if found
                {
//                    auto outit = out.find(Qout(matit)); /// see if current OUTgoing symmetry sector already is present in output matrix
//                    if (outit==out.end()) out.emplace(Qout(matit),QMat(matit).t() * ((init->second) * QMat(matit)) ); /// if not, create and store current result
//                    else outit->second += QMat(matit).t() * ((init->second) * QMat(matit)); /// if yes, add current result to already present contribution
                    auto outit = out.lower_bound(Qout(matit)); /// see if current OUTgoing symmetry sector already is present in output matrix
                    if (outit==out.end() || outit->first!=Qout(matit)) out.emplace_hint(outit,Qout(matit),QMat(matit).t() * ((init->second) * QMat(matit)) ); /// if not, create and store current result
                    else outit->second += QMat(matit).t() * ((init->second) * QMat(matit)); /// if yes, add current result to already present contribution
                }
            }
        }
    }
}

/**< This version is meant to be used in TMEigs only, where the memory of out is preallocated for speed */
template<typename KT,typename VTMPS,typename VTmat>
void
ApplyTMRight(const MPSBlockMat<KT,VTMPS>& A, const BlockDiagMat<KT,VTmat>& in, BlockDiagMat<KT,VTmat>& out)
{
    /**< ATTENTION: IF OUT IS NOT COMPLETELY EMPTY, MAKE SURE IT IS INITIALIZED TO ZERO */
    for (const auto Ait : A)
    {
        if (!Ait.empty()) /// see if there are any symmetry sectors for current physical index s
        {
            for (const auto& matit : Ait) /// loop through symmetry sectors of A_s
            {
                const auto init = in.find(Qout(matit)); /// see if current OUTgoing symmetry sector is present in input matrix
                if (init!=in.end()) /// only do something if found
                {
//                    auto outit = out.find(Qin(matit)); /// see if current INgoing symmetry sector already is present in output matrix
//                    if (outit==out.end()) out.emplace(Qin(matit), (QMat(matit) * init->second) * QMat(matit).t() ); /// if not, create and store current result
//                    else outit->second += (QMat(matit) * init->second) * QMat(matit).t(); /// if yes, add current result to already present contribution
                    auto outit = out.lower_bound(Qin(matit)); /// see if current INgoing symmetry sector already is present in output matrix
                    if (outit==out.end() || outit->first!=Qin(matit)) out.emplace_hint(outit,Qin(matit), (QMat(matit) * init->second) * QMat(matit).t() ); /// if not, create and store current result
                    else outit->second += (QMat(matit) * init->second) * QMat(matit).t(); /// if yes, add current result to already present contribution
                }
            }
        }
    }
}

/**< Mixed Transfer Matrix Operations ****************************************************************************************/

/**< Apply mixed transfer matrix, generated by A (lower) and B (upper), to the left onto unity (i.e. only \sum_s B_s' * A_s) */
template<typename KT,typename VT>
BlockDiagMat<KT,VT>
ApplyTMmixedLeft(const MPSBlockMat<KT,VT>& A, const MPSBlockMat<KT,VT>& B)
{
    assert(A.size() == B.size()); /// check if A and B have same dimensions
    BlockDiagMat<KT,VT> out;

/// TODO (valentin#1#2015-04-19): switch to iterators
    for (uint s=0; s<A.size(); ++s) /// loop through physical indices
    {
        if (!A[s].empty() && !B[s].empty()) /// check if both MPS matrices have entries for current phys. index
        {
            for (const auto Amatit : A[s]) /// for current i, loop through symmetry sectors of A_s (Amatit has type BlockMat)
            {
                const auto Bmatit = B[s].find(Qin(Amatit)); /// find corresponding symmetry sector in B_s (Bmatit has type ITERATOR of BlockMat)
                if (Bmatit!=B[s].end()) /// only do something if it is found in B_s
                {
//                    auto outit = out.find(Qout(Amatit)); /// look for current symmetry sector (result of B_s' * A_s has OUTgoing symmetry label of A_s)
//                    if (outit==out.end()) out.emplace(Qout(Amatit) , QMat(*Bmatit).t() * QMat(Amatit) );
//                    else outit->second += QMat(*Bmatit).t() * QMat(Amatit);
                    auto outit = out.lower_bound(Qout(Amatit)); /// look for current symmetry sector (result of B_s' * A_s has OUTgoing symmetry label of A_s)
                    if (outit==out.end() || outit->first!=Qout(Amatit)) out.emplace_hint(outit,Qout(Amatit) , QMat(*Bmatit).t() * QMat(Amatit) );
                    else outit->second += QMat(*Bmatit).t() * QMat(Amatit);
                }
            }
        }
    }
    return out;
}

/**< Apply mixed transfer matrix, generated by A (lower) and B (upper), to the left onto unity (i.e. only \sum_s A_s * B_s') */
template<typename KT,typename VT>
BlockDiagMat<KT,VT>
ApplyTMmixedRight(const MPSBlockMat<KT,VT>& A, const MPSBlockMat<KT,VT>& B)
{
    assert(A.size() == B.size()); /// check if A and B have same dimensions
    BlockDiagMat<KT,VT> out;

/// TODO (valentin#1#2015-04-19): switch to iterators
    for (uint s=0; s<A.size(); ++s) /// loop through physical indices
    {
        if (!A[s].empty() && !B[s].empty()) /// check if both MPS matrices have entries for current phys. index
        {
            for (const auto Amatit : A[s]) /// for current i, loop through symmetry sectors of A_s (Amatit has type BlockMat)
            {
                /// find corresponding symmetry sector in B_s (Bmatit has type ITERATOR of BlockMat),
                /// actually we should match Qout, but that's not possible. However we can exploit that if the two
                /// Qouts for the same s are equal, so are the Qins, so let's again look for these.
                const auto Bmatit = B[s].find(Qin(Amatit));
                if (Bmatit!=B[s].end()) /// only do something if it is found in B_s
                {
//                    auto outit = out.find(Qin(Amatit)); /// look for current symmetry sector (result of B_s' * A_s has OUTgoing symmetry label of A_s)
//                    if (outit==out.end()) out.emplace(Qin(Amatit) , QMat(Amatit) * QMat(*Bmatit).t());
//                    else outit->second += QMat(Amatit) * QMat(*Bmatit).t();
                    auto outit = out.lower_bound(Qin(Amatit)); /// look for current symmetry sector (result of B_s' * A_s has OUTgoing symmetry label of A_s)
                    if (outit==out.end() || outit->first!=Qin(Amatit)) out.emplace_hint(outit,Qin(Amatit) , QMat(Amatit) * QMat(*Bmatit).t());
                    else outit->second += QMat(Amatit) * QMat(*Bmatit).t();
                }
            }
        }
    }
    return out;
}
/**< DECOMPOSITIONS AND TRUNCATION =========================================================================================================== */
/// used for split2() to store info about newly generated symmetry sectors
template<typename KT,typename VT>
struct sectordata
{
    typedef typename MPSBlockMat<KT,VT>::mciter mciter; /// const_iterator through BlockMat inside a MPSBlockMat
    /// define tuple of left phys. index, right phys. index iterator to corresponding element (matrix within a BlockMat)
    /// we will need it to fill the matrix to decompose via SVD with the corresponding entries
    typedef typename std::tuple<uint,uint,mciter> ph_ind_iter; /// triple to store left phys. index, right phys. index, iterator to matrix within a BlockMat

    std::vector<ph_ind_iter> v_phys_ind; /// stores tuples of (left phys. index, right phys. index, iterator to matrix within a BlockMat)
    /// maps phys. index to pair of matrix dimension and corresponding symmetry sector
    /// we will fill two separate with all contributing left and right phys. indices and their corresponding sizes/QNs
    std::map<uint,std::pair<uint,const KT*> > sizes_left,sizes_right; /// use pointer to KeyType to avoid copying
};

template<typename KT,typename VT,uint N>
void
split(const MPSBlockMat<KT,VT>& psi, MPSBlockMat<KT,VT>& A, MPSBlockMat<KT,VT>& B, BlockLambda<KT>& lam, const ItoKey<N>& I2K)
{
    assert(psi.GetNSites()>N);
    assert(psi.GetLocalDim() == I2K.GetLocalDim());

    std::map<KT,sectordata<KT,VT> > QNmap;
    uint left,right;
    uint dimright=pow(I2K.GetLocalDim(),psi.GetNSites()-N);

    A = MPSBlockMat<KT,VT>(I2K.GetLocalDim(),N);
    B = MPSBlockMat<KT,VT>(I2K.GetLocalDim(),psi.GetNSites()-N);
    lam.clear();

    typename MPSBlockMat<KT,VT>::mciter iter;

    /// determine all possible QN (also newly generated ones) that can occur at the bond where the split is happening
    for (uint s=0; s<psi.size(); ++s)
    {
        left=s/dimright;
        right=s%dimright;

        for (iter=psi[s].begin(); iter!=psi[s].end(); ++iter)
        {
            sectordata<KT,VT>& sec(QNmap[iter->first+I2K[left]]); /// if current QN (iter->first) is not yet present in QNmap, generate corresponding sectordata object with std.c'tor

            /// for current center QN, collect all contributing elements and their sizes and (left or right) QNs
            sec.v_phys_ind.emplace_back(left,right,iter);
            sec.sizes_left.emplace(left,make_pair(QMat(*iter).n_rows,&Qin(*iter)));
            sec.sizes_right.emplace(right,make_pair(QMat(*iter).n_cols,&Qout(*iter)));
        }
    }

    /// maps of phys. indices to spans of indices in the corresponding matrices
    /// we need maps yet again, since we ASSOCIATIVELY need to access the spans via phys. indices later on
    std::map<uint,std::pair<span,const KT*> > spans_left,spans_right;
    uint ml_tot,mr_tot,lamct;
    Mat<VT> U,V;
    RVecType lamtmp;

    for (const auto& qnit : QNmap)
    {
//        cout<<"processing "<<qnit.first<<endl;
        const sectordata<KT,VT>& sec(qnit.second);
        spans_left.clear();
        spans_right.clear();
        ml_tot=0;
        mr_tot=0;

        /// construct maps of spans for left and right indices
        /// we can only do that now, since before we didn't know yet, which phys. indices actually contribute
        /// This way we can avoid decomposing matrices that are actually mostly filled with zeros (for the phys. indices that don't contribute)
        for (const auto& mlit : sec.sizes_left)
        {
            /// mlit.first = left phys. index
            /// mlit.second.first = matrix size
            /// mlit.second.second = pointer to QN
            spans_left.emplace_hint(spans_left.end(),mlit.first,std::make_pair(span(ml_tot,ml_tot+mlit.second.first-1),mlit.second.second));
            ml_tot+=mlit.second.first;
        }
        for (const auto& mrit : sec.sizes_right)
        {
            /// mrit.first = right phys. index
            /// mrit.second.first = matrix size
            /// mrit.second.second = pointer to QN
            spans_right.emplace_hint(spans_right.end(),mrit.first,std::make_pair(span(mr_tot,mr_tot+mrit.second.first-1),mrit.second.second));
            mr_tot+=mrit.second.first;
        }

        Mat<VT> M(ml_tot,mr_tot,fill::zeros);
        /// now we have all the spans for creating the matrix from the vector of contributing elements
        for (const auto& vit : sec.v_phys_ind)
        {
            left=get<0>(vit);
            right=get<1>(vit);
            /// here we need to access the spans associatively, hence the maps
            M(spans_left[left].first,spans_right[right].first) = QMat(*get<2>(vit));
        }

        /// actually perform SVD (which is the main reason for all this hustle :-) )
        svd_econ(U,lamtmp,V,M);

        /// determine if some of the Schmidt-Values are already ~zero and discard (when there has been no mixing between left and right due to e.g. a quantum gate)
        lamct=lamtmp.n_elem;
/// TODO (valentin#1#2015-05-13): Which threshold to use here?
        while (lamct>1 && lamtmp(lamct-1)<1e-14)--lamct;
        if (lamct<lamtmp.n_elem)
        {
/// TODO (valentin#1#2015-05-13): consider resize instead of recreating with head and head_cols
            lamtmp = lamtmp.head(lamct);
            U = U.head_cols(lamct);
            V = V.head_cols(lamct);
        }


        /// put results in respective left and right MPSBlockMatrices A and B and BlockLambda lam
        lam.emplace_hint(lam.end(),qnit.first,lamtmp);
        for (const auto& mlit : spans_left)
        {
            /// qnit.first = current center QN
            /// mlit.first = left phys. index
            /// mlit.first.first = span corresponding to left phys. index
            /// mlit.second.second = pointer to corresponding left QN
            if (!A[mlit.first].emplace(*mlit.second.second,QMatPair<KT,VT>(qnit.first,U(mlit.second.first,span::all))).second )
            { cerr<<"could not insert ("<<*mlit.second.second<<","<<qnit.first<<") into A["<<mlit.first<<"]"<<endl; }
        }
        for (const auto& mrit : spans_right)
        {
            /// qnit.first = current center QN
            /// mlit.first = right phys. index
            /// mlit.first.first = span corresponding to right phys. index
            /// mlit.second.second = pointer to corresponding right QN
            if (!B[mrit.first].emplace(qnit.first,QMatPair<KT,VT>(*mrit.second.second,V(mrit.second.first,span::all).t())).second)
            { cerr<<"could not insert ("<<qnit.first<<","<<*mrit.second.second<<") into B["<<mrit.first<<"]"<<endl;}
        }
    }
}


//template<typename KT,typename VT, uint N>
//void
//truncate(MPSBlockMat<KT,VT>& A, MPSBlockMat<KT,VT>& B, BlockLambda<KT>& lam, const ItoKey<N>& I2K, uint mmax, Real lamthresh=1e-10)
//{
//
//
//}

template<typename KT,typename VT, uint N>
void
truncate(MPSBlockMat<KT,VT>& A, MPSBlockMat<KT,VT>& B, BlockLambda<KT>& lam, const ItoKey<N>& I2K, uint mmax, Real lamthresh=1e-10)
{
    assert(A.GetNSites()==N);
    bool vectoobig=false,valtoosmall=false;

    typedef typename BlockLambda<KT>::iterator lam_iter_type;
    std::vector<std::pair<uint,lam_iter_type> > size_vec;
    typedef typename std::vector<std::pair<uint,lam_iter_type> >::iterator size_iter_type;

    lam_iter_type lamit;
    size_iter_type sizit;


    typedef std::pair<Real,size_iter_type> literpair;
    std::vector<literpair> lamvec;
    lamvec.reserve(lam.GetTotalSize());
//
//    A.ShowDims("A");
//    B.ShowDims("B");

/// TODO (valentin#1#2015-05-13): try not to use two loops through lam and to grow size_vec along with looping through lam
    size_vec.reserve(lam.size());
    for (lamit=lam.begin();lamit!=lam.end();++lamit) size_vec.emplace_back(0,lamit);
    for (lamit=lam.begin(),sizit = size_vec.begin();lamit!=lam.end();++lamit,++sizit)
    {
        for(const auto& vit : lamit->second)
        {
            if (vit>lamthresh) lamvec.emplace_back(vit,sizit);
            else valtoosmall=true;
        }
    }
    vectoobig = lamvec.size() > mmax;

    if (vectoobig || valtoosmall)
    {
        if (vectoobig)
        {
            std::sort(lamvec.begin(),lamvec.end(),[](const literpair& lhs,const literpair& rhs) -> bool {return lhs.first > rhs.first;});
            lamvec.resize(mmax);
        }

//        cout<<"lam"<<endl<<lam;
//        for (const auto& vit : lamvec) cout<<vit.second->second->first<<": "<<vit.first<<endl;
//        cout<<"|lam|="<<lam.GetTotalSize()<<", |lamvec|="<<lamvec.size()<<", mmax="<<mmax<<endl;

        /// determine new bond dimensions for each QN
        /// since we also already stored the iterators to the
        for (const auto& vit : lamvec)
        {
            /// vit.first = lam value
            /// vit.second = size_vec iterator
            /// vit.second->first = size
            /// vit.second->second = lam_iter_type
            ++(vit.second->first);
        }
        uint siz;
        for (const auto& vit : size_vec)
        {
            /// vit.first = size
            /// vit.second = lam_iter_type
            /// vit.second->first = KT
            /// vit.second->second = lam in sector of KT
            siz = vit.first;
            const KT& QN = vit.second->first;

//            cout<<QN<<": "<<siz<<endl;
            assert(siz <= vit.second->second.n_elem);
            if (siz < vit.second->second.n_elem) /// only do something if the size has shrunk
            {
                if (siz > 0) /// if smaller, but still finite, perform truncation
                {
//                    cout<<"truncate "<<QN<<endl;
                    vit.second->second.resize(siz);
                    for (uint s=0;s<A.size();++s)
                    {
                        auto it = A[s].find(QN - I2K[s]);
                        if (it!=A[s].end()) QMat(*it).resize(QMat(*it).n_rows,siz);
//                        if (it!=A[s].end()) {QMat(*it).resize(QMat(*it).n_rows,siz);cout<<"resizing ("<<Qin(*it)<<","<<Qout(*it)<<") of A["<<s<<"]"<<endl;}
//                        else cout<<QN - I2K[s]<<" not found in A["<<s<<"]"<<endl;
                    }
                    for (uint s=0;s<B.size();++s)
                    {
                        auto it = B[s].find(QN);
                        if (it!=B[s].end()) QMat(*it).resize(siz,QMat(*it).n_cols);
//                        if (it!=B[s].end()) {QMat(*it).resize(siz,QMat(*it).n_cols);cout<<"resizing ("<<Qin(*it)<<","<<Qout(*it)<<") of B["<<s<<"]"<<endl;}
//                        else cout<<QN<<" not found in A["<<s<<"]"<<endl;
                    }
                }
                else /// if completely gone, erase from containers
                {
//                    cout<<"erase "<<QN<<endl;
//                    lam.erase(vit.second);

                    for (uint s=0;s<A.size();++s)
                    {
                        A[s].erase(QN - I2K[s]);
//                        if (A[s].erase(QN - I2K[s])>0 ) cout<<"erased ("<<QN - I2K[s]<<","<<QN<<") from A["<<s<<"]"<<endl;
//                        else cout<<QN - I2K[s]<<" not found in A["<<s<<"]"<<endl;
                    }
                    for (uint s=0;s<B.size();s++)
                    {
                        B[s].erase(QN);
//                        if(B[s].erase(QN)>0) cout<<"erased ("<<QN<<","<<QN + I2K[s]<<") from B["<<s<<"]"<<endl;
//                        else cout<<QN<<" not found in B["<<s<<"]"<<endl;
                    }

                    if(lam.erase(vit.second)==vit.second)cerr<<QN<<" should have been erased from lam!"<<endl;
//                    else cout<<"erased "<<QN<<" from lam"<<endl;

                }
            }
        }
    }
    else cout<<"nothing to do"<<endl;
}


/**< ===================================================================================================================================== */
/**< random uMPS */
template<typename VT,typename KT,uint N>
MPSBlockMat<KT,VT>
RandUMPS(const map<KT,uint>& dimmap,const ItoKey<N>& I2K)
{
    srand(time(NULL));
    MPSBlockMat<KT,VT> out(I2K.GetLocalDim(),N);
    for (uint s=0; s<I2K.size(); ++s)
    {
        for (const auto& lkit : dimmap)
        {
            auto rkit = dimmap.find(lkit.first + I2K[s]);
            if (rkit!=dimmap.end()) out[s].emplace(lkit.first, QMatPair<KT,VT>(rkit->first, randn<Mat<VT> >(lkit.second,rkit->second)) );
        }
    }
    return out;
}

template<typename VT,typename KT>
void
CheckOrtho(const MPSBlockMat<KT,VT>& A, const BlockLambda<KT>& lam, dirtype dir)
{
    Real chk1=0,chk2=0;
    BlockDiagMat<KT,VT> L,R,lammat;
    std::string str;
/// TODO (valentin#1#2015-05-11): Implement operator+(matrix,vector) to subtract a diagonal matrix from a normal matrix
    switch (dir)
    {
    case l:
        str = "left";
        lammat = pow(lam,2);
        L = ApplyTMLeft(A);
        R = ApplyTMRight(A<<lam) - lammat;
        for (const auto& it : L) chk1 += pow(norm(it.second - eye(it.second.n_rows,it.second.n_cols),"fro"),2)/it.second.n_elem;
        for (const auto& it : R) chk2 += pow(norm(it.second,"fro"),2)/it.second.n_elem;
//        for (const auto& it : L) chk1 += dot(it.second - eye(it.second.n_rows,it.second.n_cols))/it.second.n_elem;
//        for (const auto& it : R) chk2 += dot(it.second)/it.second.n_elem;

        break;
    case r:
        str = "right";
        lammat = pow(lam,2);
        L = ApplyTMLeft(lam>>A) - lammat;
        R = ApplyTMRight(A);
        for (const auto& it : L) chk1 += pow(norm(it.second,"fro"),2)/it.second.n_elem;
        for (const auto& it : R) chk2 += pow(norm(Mat<VT>(it.second - eye(it.second.n_rows,it.second.n_cols)),"fro"),2)/it.second.n_elem;
//        for (const auto& it : L) chk1 += dot(it.second)/it.second.n_elem;
//        for (const auto& it : R) chk2 += dot(Mat<VT>(it.second - eye(it.second.n_rows,it.second.n_cols)))/it.second.n_elem;
        break;
    case s:
        str="symmetric";
        lammat = lam;
        L = ApplyTMmixedLeft(lam>>A,A) - lammat;
        R = ApplyTMmixedRight(A<<lam,A) - lammat;
        for (const auto& it : L) chk1 += pow(norm(it.second,"fro"),2)/it.second.n_elem;
        for (const auto& it : R) chk2 += pow(norm(it.second,"fro"),2)/it.second.n_elem;
//        for (const auto& it : L) chk1 += dot(it.second)/it.second.n_elem;
//        for (const auto& it : R) chk2 += dot(it.second)/it.second.n_elem;
        break;
    case c:
        str="canonical";
        L = ApplyTMLeft(lam>>A);
        R = ApplyTMRight(A<<lam);
        for (const auto& it : L) chk1 += pow(norm(Mat<VT>(it.second - eye(it.second.n_rows,it.second.n_cols)),"fro"),2)/it.second.n_elem;
        for (const auto& it : R) chk2 += pow(norm(Mat<VT>(it.second - eye(it.second.n_rows,it.second.n_cols)),"fro"),2)/it.second.n_elem;
//        for (const auto& it : L) chk1 += dot(Mat<VT>(it.second - eye(it.second.n_rows,it.second.n_cols)))/it.second.n_elem;
//        for (const auto& it : R) chk2 += dot(Mat<VT>(it.second - eye(it.second.n_rows,it.second.n_cols)))/it.second.n_elem;
        break;
    default:
        cerr<<"wrong direction specified"<<endl;
    }
    cout<<"check "<<str<<" gauge:"<<endl;
    cout<<"left: "<<sqrt(chk1)<<", right:"<<sqrt(chk2)<<endl;
}


#endif // MPS_BLOCK_UTIL_H
