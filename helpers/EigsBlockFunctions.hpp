#ifndef EIGS_BLOCK_FUN_H_
#define EIGS_BLOCK_FUN_H_

#include <string>
//#include <typeinfo>

#include "arma_typedefs.h"
#include "MPSBlockMat.hpp"
#include "eigs.h"

using std::map;
using std::cout;
using std::endl;
using std::string;

/// Multiplication routine for eigs, which works just with C-style arrays of double.
/// BlockDiagMats therefore get vectorized by successively vectorizing all matrices for every symmetry sector
/// and concatenating them consecutively. To avoid copying as much as possible, Arma matrices get binded to the
/// respective memory sectors within these arrays to enable high level routines with quantum numbers for applying the TM
template<typename KT,typename VT>
void
TMMult(VT* invec, VT* outvec, const std::function<void (const BlockDiagMat<KT,VT>&, BlockDiagMat<KT,VT>&)>& TMfun, const map<KT,uint>& dimmap, uint D_tot)
{
    BlockDiagMat<KT,VT> in,out;
    memset(outvec,0,D_tot*(sizeof(*outvec)));

    uint pos=0;
    for (const auto& dimit : dimmap)
    {
        in.emplace_hint(in.end(),dimit.first,Mat<VT>(&invec[pos], dimit.second, dimit.second, false, true));
        out.emplace_hint(out.end(),dimit.first,Mat<VT>(&outvec[pos], dimit.second, dimit.second, false, true));
        pos += dimit.second*dimit.second;
    }
    TMfun(in,out);
}

template<typename KT>
inline Real
TMDominantEig(const MPSBlockMat<KT,Real>& MPS, BlockDiagMat<KT,Real>& V, dirtype dir, Real tol=1e-14, int maxit=500, string mode="LR")
{
    map<KT,uint> dimmap = MPS.GetAllSizes(); /// get all possible symmetry sectors present as ingoing or outgoing from MPS and their dimensions
    uint D_tot=0;
    for (const auto& it : dimmap) D_tot += it.second*it.second; /// sum up all dimensions to get maximum dimension (for linear vector)
    return TMDominantEig(MPS,dimmap,D_tot,V,dir,tol,maxit,mode);
}

template<typename KT>
Real
TMDominantEig(const MPSBlockMat<KT,Real>& MPS, const map<KT,uint>& dimmap, uint D_tot, BlockDiagMat<KT,Real>& V, dirtype dir, Real tol=1e-14, int maxit=500, string mode="LR")
{
    if (mode!="LR" && mode!="LM")
    {
        cerr<<"wrong mode "<<mode<<endl;
        abort();
    }
    /// function handle for the actual routine for applying the TM frm the left/right onto some BlockDiagMat
    std::function<void (const BlockDiagMat<KT,Real>&, BlockDiagMat<KT,Real>&)> TMfun;
    if (dir==l) TMfun=[&MPS](const BlockDiagMat<KT,Real>& in, BlockDiagMat<KT,Real>& out) -> void {ApplyTMLeft(MPS,in,out);};
    else if (dir==r) TMfun=[&MPS](const BlockDiagMat<KT,Real>& in, BlockDiagMat<KT,Real>& out) -> void {ApplyTMRight(MPS,in,out);};
    else
    {
        cerr<<"wrong direction specified"<<endl;
        abort();
    }

    /// actual calculation of the dominant eigenpair of the TM
    CVecType valtmp;
    CMatType Vtmp;
    eigs_rn([&dimmap,D_tot,&TMfun](Real* invec, Real* outvec)->void {TMMult(invec,outvec,TMfun,dimmap,D_tot);},D_tot,valtmp,Vtmp,1,mode,tol,maxit);

    /// analyze and postedit the dominant eigenpair
    V.clear();
    if (imag(valtmp(0))>2*tol) cerr<<"Warning: dominant eigenvalue is complex: "<<valtmp(0)<<endl;
    if (norm(imag(Vtmp.col(0)))>2*tol*D_tot) cerr<<"Warning: dominant eigenvector is complex"<<endl;
    RVecType Vvec = real(Vtmp.col(0));
    Real val = real(valtmp(0));

    uint pos=0;
//    Real tr=0; /// compute overall trace to make eigenmat. positive in the end
    for (const auto& it : dimmap)
    {
        RMatType tmp(&Vvec.memptr()[pos],it.second,it.second);
//        tmp = tmp + tmp.t(); /// hermiticize (apparently not necessary when using quantum numbers!?)
        V.emplace_hint(V.end(),it.first,tmp);
//        V.emplace_hint(V.end(),it.first,RMatType(&Vvec.memptr()[pos],it.second,it.second));
        pos += it.second*it.second;
//        tr += trace(tmp);
    }
    V/=trace(V);
//    for (auto& it : V) it.second/=tr; /// make positive

    if (dir==l) DOUT("left: "<<abs(val)<<", "<<val<<": "<<norm(ApplyTMLeft(MPS,V) - val*V));
    else if (dir==r) DOUT("right: "<<abs(val)<<", "<<val<<": "<<norm(ApplyTMRight(MPS,V) - val*V));

    return val;
}

template<typename KT>
inline CVecType
TMEigs(const MPSBlockMat<KT,Real>& MPS, std::vector<BlockDiagMat<KT,Complex> >& V, dirtype dir, uint nev, Real tol=1e-14, int maxit=500, string mode="LR")
{
    map<KT,uint> dimmap = MPS.GetAllSizes(); /// get all possible symmetry sectors present as ingoing or outgoing from MPS and their dimensions
    uint D_tot=0;
    for (const auto& it : dimmap) D_tot += it.second*it.second; /// sum up all dimensions to get maximum dimension (for linear vector)
    return TMEigs(MPS,dimmap,D_tot,V,dir,nev,tol,maxit,mode);
}



template<typename KT>
CVecType
TMEigs(const MPSBlockMat<KT,Real>& MPS, const map<KT,uint>& dimmap, uint D_tot, std::vector<BlockDiagMat<KT,Complex> >& V, dirtype dir, uint nev, Real tol=1e-14, int maxit=500, string mode="LR")
{
    if (mode!="LR" && mode!="LM")
    {
        cerr<<"wrong mode "<<mode<<endl;
        abort();
    }
    /// function handle for the actual routine for applying the TM frm the left/right onto some BlockDiagMat
    std::function<void (const BlockDiagMat<KT,Real>&, BlockDiagMat<KT,Real>&)> TMfun;
    if (dir==l) TMfun=[&MPS](const BlockDiagMat<KT,Real>& in, BlockDiagMat<KT,Real>& out) -> void {ApplyTMLeft(MPS,in,out);};
    else if (dir==r) TMfun=[&MPS](const BlockDiagMat<KT,Real>& in, BlockDiagMat<KT,Real>& out) -> void {ApplyTMRight(MPS,in,out);};
    else
    {
        cerr<<"wrong direction specified"<<endl;
        abort();
    }


    /// actual calculation of the dominant eigenpair of the TM
    CVecType vals;
    CMatType Vtmp;
//    cout<<"mode="<<mode<<endl;
//    eigs_rn(MultOPx,D_tot,vals,Vtmp,nev,mode,tol,maxit);
    eigs_rn([&dimmap,D_tot,&TMfun](Real* invec, Real* outvec)->void {TMMult(invec,outvec,TMfun,dimmap,D_tot);},D_tot,vals,Vtmp,nev,mode,tol,maxit);

    nev = vals.size();
    /// analyze and postedit the dominant eigenpair
    V.clear();
    V.resize(nev);

    /// reorder
    if (nev>1)
    {
        uvec order;
        if (mode=="LR") order = sort_index(real(vals),"descend");
        else if (mode=="LM") order = sort_index(abs(vals),"descend");
        vals = vals(order);
        Vtmp = Vtmp.cols(order);
    }

    for (uint n=0; n<nev; ++n)
    {
        uint pos=0;
        for (const auto& it : dimmap)
        {
            CMatType tmp(&(Vtmp.colptr(n))[pos],it.second,it.second);
            V[n].emplace_hint(V[n].end(),it.first,tmp);
            pos += it.second*it.second;
        }
        if (dir==l) DOUT("left: "<<abs(vals(n))<<", "<<vals(n)<<": "<<norm(ApplyTMLeft(MPS,V[n]) - vals(n)*V[n]));
        else if (dir==r) DOUT("right: "<<abs(vals(n))<<", "<<vals(n)<<": "<<norm(ApplyTMRight(MPS,V[n]) - vals(n)*V[n]));
    }
    DOUT("");
    return vals;
}

template<typename KT,typename VT>
MPSBlockMat<KT,VT>
OrthoUMPS(BlockLambda<KT>& lam, const MPSBlockMat<KT,VT>& A, dirtype dir, double tol=1e-14, int maxit=500)
{
    map<KT,uint> dimmap = A.GetAllSizes(); /// get all possible symmetry sectors present as ingoing or outgoing from MPS and their dimensions
    uint D_tot=0;
    for (const auto& it : dimmap) D_tot += it.second*it.second; /// sum up all dimensions to get maximum dimension (for linear vector)

    BlockDiagMat<KT,VT> Vl,Vr;
    Real vall = TMDominantEig(A,dimmap,D_tot,Vl,l,tol,maxit);
    Real valr = TMDominantEig(A,dimmap,D_tot,Vr,r,tol,maxit);

    if (abs(vall-valr)>1e-10) cerr<<"valr and valr differ by "<<abs(vall-valr)<<endl;
//    DOUT("left ("<<vall<<"): "<<norm(ApplyTMLeft(A,Vl) - vall*Vl)<<", hermitian: "<<norm(Vl - trans(Vl)));
//    DOUT("right ("<<valr<<"): "<<norm(ApplyTMRight(A,Vr) - valr*Vr)<<", hermitian: "<<norm(Vr - trans(Vr)));

    Real nrm = 0.5*(vall + valr);
    BlockDiagMat<KT,Real> matl,matr;
    typename BlockDiagMat<KT,VT>::const_iterator vlit,vrit;

    Mat<VT> Ul,Ur,U,V,X,Y,ltmp,rtmp;
    RVecType Dl,Dr,lamtmp,Dlsq,Drsq,Dlsqi,Drsqi;
    assert(Vl.size()==Vr.size());
    Real eigthresh = 1e-10;

    Real nrml=0,nrmr=0,lamnrm2=0;
    lam.clear();
    assert(Vl.size() == Vr.size());
//    cout<<"tr(Vl): "<<trace(Vl)<<endl;
//    cout<<"tr(Vr): "<<trace(Vr)<<endl;
//    for (vlit = Vl.begin(),vrit = Vr.begin(); vlit!=Vl.end(); vlit++,vrit++)
//    {
//        eig_sym(Dl,Ul,vlit->second);
//        eig_sym(Dr,Ur,vrit->second);
//        Dl.print("Dl");
//        Dr.print("Dr");
//    }
    for (vlit = Vl.begin(),vrit = Vr.begin(); vlit!=Vl.end(); vlit++,vrit++)
    {
        nrml=norm(vlit->second,"fro");
        nrmr=norm(vrit->second,"fro");
        assert(vlit->first == vrit->first);
        if (nrml>eigthresh && nrmr>eigthresh)
        {
//            DOUT(vlit->first);
            eig_sym(Dl,Ul,vlit->second);
            eig_sym(Dr,Ur,vrit->second);
//            Dl.print("Dl");
//            Dr.print("Dr");
            uword indl=0,indr=0,nsiz;
            while (indl<Dl.n_elem && Dl(indl)<eigthresh)++indl;
            while (indr<Dr.n_elem && Dr(indr)<eigthresh)++indr;

            if (indl>0)
            {
                nsiz = Dl.n_elem-indl;
                if (nsiz<1) {cerr<<"Dl: all eigenvalues zero"<<endl;Dl.print("Dl");}
                Dl = Dl.tail(nsiz);
                Ul = Ul.tail_cols(nsiz);
            }
            if (indr>0)
            {
                nsiz = Dr.n_elem-indr;
                if (nsiz<1) {cerr<<"Dr: all eigenvalues zero"<<endl;Dr.print("Dr");}
                Dr = Dr.tail(nsiz);
                Ur = Ur.tail_cols(nsiz);
            }
            Dlsq = sqrt(Dl);
            Drsq = sqrt(Dr);

            X = Ur<<Drsq;
            Y = Dlsq>>Mat<VT>(Ul.t());
//            DOUT("decomposition of Vr: "<<norm(vrit->second - X*X.t(),"fro"));
//            DOUT("decomposition of Vl: "<<norm(vlit->second - Y.t()*Y,"fro"));

            svd(U,lamtmp,V,Y*X);
            lamnrm2 += dot(lamtmp,lamtmp);
            lam.emplace_hint(lam.end(),vlit->first,lamtmp);

            matl.emplace_hint(matl.end(),vlit->first,(Mat<VT>(V.t()) >> Drsq) * Ur.t());
            matr.emplace_hint(matr.end(),vlit->first,Ul * (Dlsq << U));

        }
    }
    lam/=sqrt(lamnrm2);

    MPSBlockMat<KT,VT> B(A.GetLocalDim(),A.GetNSites());
    BlockLambda<KT> lamsq;

    switch (dir)
    {
    case l:
        B = ((sqrt(lamnrm2/nrm) * lam) >> matl) * (A * matr);
        break;
    case r:
        B =  (matl * A) * (matr << (lam * sqrt(lamnrm2/nrm)));
        break;
    case s:
        lamsq = sqrt(lam);
        B =  ((sqrt(lamnrm2/nrm) * lamsq) >> matl) * (A* (matr << lamsq ) ) ;
        break;
    case c:
        B = (sqrt(lamnrm2/nrm) * matl) * (A* matr) ;
        break;
    default:
        cerr<<"wrong direction specified"<<endl;
    }
    return B;
}
#endif // EIGS_BLOCK_FUN_H_
