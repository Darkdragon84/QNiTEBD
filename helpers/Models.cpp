//#warning "including models from implementation"
#include "Models.h"
#ifdef _USE_SYMMETRIES_
#include "KeyTypes.hpp"
#endif // _USE_SYMMETRIES_

using namespace std;

/**< AUXILIARY GLOBAL FUNCTIONS ********************************************************************/
const std::vector<size_t> num2ditvec(size_t x, size_t d, size_t N)
{
    assert(x<std::pow(d,N));
    std::vector<size_t> v;
    v.reserve(N);
    size_t div;
    for (size_t i=N-1; i>0; --i)
    {
        div=std::pow(d,i);
        v.push_back(x/div);
        x=x%div;
    }
    v.push_back(x);
    return v;
}

size_t ditvec2num(const std::vector<size_t>& vec,size_t d)
{
    size_t x=0,i=0;
    std::vector<size_t>::const_reverse_iterator it;
    for (it=vec.rbegin(); it!=vec.rend(); ++it)
    {
        x+=(*it)*std::pow(d,i++);
    }
    return x;
}

/**< GENERAL MODEL TYPE ***************************************************************************/
const SparseOperator<double> ModelType::ImagTimeEvoOp(Real dt) const
{
    RMatType V;
    RVecType E;
    eig_sym(E,V,ham_.Full());

    RMatType tmp(V);
    for (uint i=0; i<E.n_elem; ++i) tmp.col(i)*=std::exp(-dt*E(i));
//    RMatType U = ;
    SparseOperator<double> U(tmp*RMatType(V.t()),localdim_);
    return U;
}

/**< XYZ BASE MODEL TYPE ***************************************************************************/

/// XYZ Base Model C'tor
XYZModelType::XYZModelType(double Jx, double Jy, double Jz, double hz, size_t dim, size_t Nsyms):
    ModelType(dim,Nsyms),Jx_(Jx),Jy_(Jy),Jz_(Jz),hz_(hz)
{
    if(dim!=2 && dim!=3)cerr<<"local dim "<<dim<<" not implemented"<<endl;
//    this->localdim_=dim;
    this->spin_[2]="spin-1/2 ";
    this->spin_[3]="spin-1 ";
    MakeOps();
}

void XYZModelType::ShowParams() const
{

    cout<<"localdim: "<<this->localdim_<<endl;
    cout<<"model parameters:"<<endl;
    cout<<"Ising interaction Jx="<<Jx_<<", Jy="<<Jy_<<", Jz="<<Jz_<<endl;
    cout<<"Zeeman magnetic field hz="<<hz_<<endl;
}

void XYZModelType::MakeOps()
{
    double fac=1./std::sqrt(2.);
    id_.SetDims(1,localdim_);
    for (uint i=0; i<localdim_; ++i)id_(OpKey(i,i))=1.;

    switch (localdim_)
    {
    case 2:
        sx_.SetDims(1,2);
        sx_(OpKey(0,1))=0.5;
        sx_(OpKey(1,0))=0.5;

        syi_.SetDims(1,2);
        syi_(OpKey(0,1))=-0.5;
        syi_(OpKey(1,0))=0.5;

        sz_.SetDims(1,2);
        sz_(OpKey(0,0))=0.5;
        sz_(OpKey(1,1))=-0.5;
        break;
    case 3:
        sx_.SetDims(1,3);
        sx_(OpKey(0,1))=fac;
        sx_(OpKey(1,0))=fac;
        sx_(OpKey(1,2))=fac;
        sx_(OpKey(2,1))=fac;

        syi_.SetDims(1,3);
        syi_(OpKey(0,1))=-fac;
        syi_(OpKey(1,0))=-fac;
        syi_(OpKey(1,2))=fac;
        syi_(OpKey(2,1))=fac;

        sz_.SetDims(1,3);
        sz_(OpKey(0,0))=1.;
        sz_(OpKey(2,2))=-1.;
        break;
    default:
        cerr<<"MakeOps(): dimension not implemented"<<endl;
    }
    ham_ = -Jx_*kron(sx_,sx_) + Jy_*kron(syi_,syi_) - Jz_*kron(sz_,sz_) - 0.5*hz_*(kron(sz_,id_) + kron(id_,sz_));
//    sx_.ShowElements("sx");
//    syi_.ShowElements("syi");
//    sz_.ShowElements("sz");
}

/**< XYZ MODEL TYPE ***************************************************************************/


/**< XXZ MODEL TYPE ***************************************************************************/



/**< TFI MODEL TYPE ***************************************************************************/

/**< FERMI HUBBARD MODEL TYPE *******************************************************************/

FHubModel::FHubModel(double t,double U,double mu,double hz):
    ModelType(4,2),t_(t),U_(U),mu_(mu),hz_(hz)
{
    MakeOps();
#ifdef _USE_SYMMETRIES_
    InitSyms();
#endif // _USE_SYMMETRIES_
}

void FHubModel::ShowParams () const
{
    cout<<"Hubbard type model:"<<endl;
    cout<<"localdim: "<<this->localdim_<<endl;
    cout<<"model parameters:"<<endl;
    cout<<"hopping t="<<t_<<endl;
    cout<<"onsite interaction U="<<U_<<endl;
    cout<<"chem. potential mu="<<mu_<<endl;
    cout<<"Zeeman magnetic field hz="<<hz_<<endl;
}

void FHubModel::MakeOps()
{
    id_.SetDims(1,localdim_);
    for (uint i=0; i<localdim_; ++i)id_(OpKey(i,i))=1.;

    cup_.SetDims(1,4);
    cup_(OpKey(0,2))=1.;
    cup_(OpKey(1,3))=1.;

    cdo_.SetDims(1,4);
    cdo_(OpKey(0,1))=1.;
    cdo_(OpKey(2,3))=-1.;

    nup_.SetDims(1,4);
    nup_(OpKey(2,2))=1.;
    nup_(OpKey(3,3))=1.;

    ndo_.SetDims(1,4);
    ndo_(OpKey(1,1))=1.;
    ndo_(OpKey(3,3))=1.;

    SparseOperator<double> nup2(nup_-0.5*id_), ndo2(ndo_-0.5*id_);
    SparseOperator<double> nund(nup2*ndo2),nupnd(nup_+ndo_),sz(0.5*(nup_-ndo_));
    SparseOperator<double> F(1,4);
    F(OpKey(0,0))=1.;
    F(OpKey(1,1))=-1.;
    F(OpKey(2,2))=-1.;
    F(OpKey(3,3))=1.;

//    SparseOperator<double> cdoD(herm(cdo_));
//    cdo_.ShowElements("cdo");
//    cdoD.ShowElements("cdo'");

    ham_ = -t_*(kron(herm(cup_)*F,cup_) - kron(cup_*F,herm(cup_)) + kron(herm(cdo_)*F,cdo_) - kron(cdo_*F,herm(cdo_)))
           + 0.5*U_*(kron(nund,id_) + kron(id_,nund))
           - 0.5*mu_*(kron(nupnd,id_) + kron(id_,nupnd))
           - 0.5*hz_*(kron(sz,id_) + kron(id_,sz));
}



/**< EXTENDED FERMI HUBBARD MODEL TYPE **********************************************************/
ExtFHubModel::ExtFHubModel(double t,double U, double mu,double hz, double V, double a):
    FHubModel(t,U,mu,hz),V_(V),a_(a)
{
    FHubModel::MakeOps();
    cerr<<"ATTENTION: EXTENDED MODEL OPERATORS NOT ADDED YET!"<<endl;
#ifdef _USE_SYMMETRIES_
    InitSyms();
#endif // _USE_SYMMETRIES_
}

void ExtFHubModel::ShowParams () const
{
    cout<<"Extended ";
    FHubModel::ShowParams();
    cout<<"nearest neighbor interaction V="<<V_<<endl;
    cout<<"spin-orbit coupling alpha="<<a_<<endl;
}


/**< TEST 2*Z2 MODEL TYPE **********************************/
Z2TestModel::Z2TestModel():
    ModelType(4,2)
    {
#ifdef _USE_SYMMETRIES_
        InitSyms();
        #endif // _USE_SYMMETRIES_
    };

/**< SYMMETRIES RELATED METHODS ******************************************************************/
#ifdef _USE_SYMMETRIES_

/**< XYZ MODEL TYPE ***************************************************************************/
void XYZModel::InitSyms()
{
    /// XYZ has Z2 symmetry
    if(localdim_==2)
    {
        keys_.reserve(2);
        keys_.push_back({0}); /**< spin down (state 0) has quantum number 0 */
        keys_.push_back({1}); /**< spin up (state 1) has quantum number 1 */
    }
    else cerr<<"InitSyms(): local dimension "<<localdim_<<" not implemented"<<endl;
}

/**< For Z2, group addition and subtraction are the same! */
IKey XYZModel::Plus(const IKey& lhs, const IKey& rhs) const { return (lhs.Plus(rhs)).Modulus(2);}
IKey XYZModel::Minus(const IKey& lhs, const IKey& rhs) const { return Plus(lhs,rhs);} /// for


//MPSBlockMat<IKey,Real>
//XYZModel::RandRMPSMat(size_t ml0, size_t mr0, size_t Nsecs) const//, size_t Nsites=1)
//{
//    size_t ml=ml0/2,mr=mr0/2,Nsites=1;
//    double nrm = pow(ml0*mr0,0.25);
//    srand(time(NULL));
//
//    MPSBlockMat<IKey,Real> out(this->localdim_,Nsites);
//    IKey K1({0}),K2({1});
//
//    out[0][IMPSKey(K1,this->Plus(K1,this->keys_[0]))] = randn<RMatType>(ml,mr)/nrm;
//    out[0][IMPSKey(K2,this->Plus(K2,this->keys_[0]))] = randn<RMatType>(ml,mr)/nrm;
//
//    out[1][IMPSKey(K1,this->Plus(K1,this->keys_[1]))] = randn<RMatType>(ml,mr)/nrm;
//    out[1][IMPSKey(K2,this->Plus(K2,this->keys_[1]))] = randn<RMatType>(ml,mr)/nrm;
//
//    return out;
//}

//MPSBlockMat<IKey,Complex>
//XYZModel::RandCMPSMat(size_t ml0, size_t mr0, size_t Nsecs) const//, size_t Nsites=1)
//{
//    size_t ml=ml0/2,mr=mr0/2,Nsites=1;
//    srand(time(NULL));
//
//    MPSBlockMat<IKey,Complex> out(this->localdim_,Nsites);
//    IKey K1({0}),K2({1});
//
//    out[0][IMPSKey(K1,this->Plus(K1,this->keys_[0]))] = randn<CMatType>(ml,mr);
//    out[0][IMPSKey(K2,this->Plus(K2,this->keys_[0]))] = randn<CMatType>(ml,mr);
//
//    out[1][IMPSKey(K1,this->Plus(K1,this->keys_[1]))] = randn<CMatType>(ml,mr);
//    out[1][IMPSKey(K2,this->Plus(K2,this->keys_[1]))] = randn<CMatType>(ml,mr);
//
//    return out;
//}

/**< XXZ MODEL TYPE ***************************************************************************/
void
XXZModel::InitSyms()
{
    /// XXZ has U1 symmetry (magnetization in z)
    if(localdim_==2)
    {
        keys_.reserve(2);
        keys_.push_back({-1}); /**< spin down (state 0) has quantum number -1 */
        keys_.push_back({1}); /**< spin up (state 1) has quantum number 1 */
    }
    else cerr<<"InitSyms(): dimension not implemented"<<endl;
}

IKey XXZModel::Plus(const IKey& lhs, const IKey& rhs) const { return (lhs.Plus(rhs));}
IKey XXZModel::Minus(const IKey& lhs, const IKey& rhs) const { return (lhs.Minus(rhs));}


//MPSBlockMat<IKey,Real>
//XXZModel::RandRMPSMat(size_t ml0, size_t mr0, size_t Nsecs) const//, size_t Nsites=1)
//{
//    ivec secs(Nsecs+1); /// create one more symmetry sector, s.t. |0> gets 2:N and |1> gets 1:(N-1)
//
//    if ( Nsecs%2 == 0) secs=linspace<ivec>(-Nsecs/2,Nsecs/2,Nsecs+1);
//    else secs = linspace<ivec>(-(int)(Nsecs/2+1),Nsecs/2,Nsecs+1);
//
//    /// uniform distribution of matrix dimensions into symmetry sectors
//    ivec ml = (ml0/(Nsecs + 1))*ones<ivec>(Nsecs+1);
//    ivec mr = (mr0/(Nsecs + 1))*ones<ivec>(Nsecs+1);
//
////    /// normal distribution of matrix dimensions into symmetry sectors
////    double s = (double)Nsecs/4.;
////    double nrm = sqrt(2*datum::pi)*s*std::erf(-secs(0)/(sqrt(2)*s));
////    vec expvec = exp(-pow(conv_to<vec>::from(secs),2)/(2*s*s))/nrm;
////
////    ivec ml = conv_to<ivec>::from(round(ml0*expvec));
////    ivec mr = conv_to<ivec>::from(round(mr0*expvec));
//
//    double nrm=sqrt(2*sqrt(mean(ml)*mean(mr)));
////    double nrm=0;
////    for (uint i=0;i<Nsecs;++i) nrm +=sqrt(ml(i)*mr(i)) + sqrt(ml(i+1)*mr(i+1));
////    nrm = sqrt(nrm);
//
//    MPSBlockMat<IKey,Real> out(this->localdim_,1);
//    IKey K;
//    IMPSKey MK;
//    for (uint i=0; i<Nsecs; ++i)
//    {
//        K = IKey({secs(i+1)});
//        MK = IMPSKey(K,this->Plus(K,this->keys_[0]));
//        out[0][MK] = randn(ml(i+1),mr(i))/nrm;
//
//        K = IKey({secs(i)});
//        MK = IMPSKey(K,this->Plus(K,this->keys_[1]));
//        out[1][MK] = randn(ml(i),mr(i+1))/nrm;
//    }
//
//    return out;
//}
//
//MPSBlockMat<IKey,Complex>
//XXZModel::RandCMPSMat(size_t ml0, size_t mr0, size_t Nsecs) const//, size_t Nsites=1)
//{
//    ivec secs(Nsecs+1); /// create one more symmetry sector, s.t. |0> gets 2:N and |1> gets 1:(N-1)
//
//    if ( Nsecs%2 == 0) secs=linspace<ivec>(-Nsecs/2,Nsecs/2,Nsecs+1);
//    else secs = linspace<ivec>(-(int)(Nsecs/2+1),Nsecs/2,Nsecs+1);
//
//    /// uniform distribution of matrix dimensions into symmetry sectors
//    ivec ml = (ml0/(Nsecs + 1))*ones<ivec>(Nsecs+1);
//    ivec mr = (mr0/(Nsecs + 1))*ones<ivec>(Nsecs+1);
//
////    /// normal distribution of matrix dimensions into symmetry sectors
////    double s = (double)Nsecs/4.;
////    double nrm = sqrt(2*datum::pi)*s*std::erf(-secs(0)/(sqrt(2)*s));
////    vec expvec = exp(-pow(conv_to<vec>::from(secs),2)/(2*s*s))/nrm;
////
////    ivec ml = conv_to<ivec>::from(round(ml0*expvec));
////    ivec mr = conv_to<ivec>::from(round(mr0*expvec));
//
//    double nrm=sqrt(2*sqrt(mean(ml)*mean(mr)));
//
//    MPSBlockMat<IKey,Complex> out(this->localdim_,1);
//    IKey K;
//    IMPSKey MK;
//    for (uint i=0; i<Nsecs; ++i)
//    {
//        K = IKey({secs(i+1)});
//        MK = IMPSKey(K,this->Plus(K,this->keys_[0]));
//        out[0][MK] = randn<Mat<Complex> >(ml(i+1),mr(i))/nrm;
//
//        K = IKey({secs(i)});
//        MK = IMPSKey(K,this->Plus(K,this->keys_[1]));
//        out[1][MK] = randn<Mat<Complex> >(ml(i),mr(i+1))/nrm;
//    }
//
//    return out;
//}

/**< TFI MODEL TYPE ***************************************************************************/
void
TransverseIsingModel::InitSyms()
{
    /// TFI has Z2 symmetry
    keys_.reserve(2);
    keys_.push_back({0}); /**< spin down (state 0) has quantum number 0 */
    keys_.push_back({1}); /**< spin up (state 1) has quantum number 1 */
}

IKey TransverseIsingModel::Plus(const IKey& lhs, const IKey& rhs) const { return (lhs.Plus(rhs)).Modulus(2);}
IKey TransverseIsingModel::Minus(const IKey& lhs, const IKey& rhs) const { return TransverseIsingModel::Plus(lhs,rhs);}

//MPSBlockMat<IKey,Real>
//TransverseIsingModel::RandRMPSMat(size_t ml0, size_t mr0, size_t Nsecs) const//, size_t Nsites=1)
//{
//    size_t ml=ml0/2,mr=mr0/2,Nsites=1;
//    double nrm = pow(ml0*mr0,0.25);
//    srand(time(NULL));
//
//    MPSBlockMat<IKey,Real> out(this->localdim_,Nsites);
//    IKey K1({0}),K2({1});
//
//    out[0][IMPSKey(K1,this->Plus(K1,this->keys_[0]))] = randn<RMatType>(ml,mr)/nrm;
//    out[0][IMPSKey(K2,this->Plus(K2,this->keys_[0]))] = randn<RMatType>(ml,mr)/nrm;
//
//    out[1][IMPSKey(K1,this->Plus(K1,this->keys_[1]))] = randn<RMatType>(ml,mr)/nrm;
//    out[1][IMPSKey(K2,this->Plus(K2,this->keys_[1]))] = randn<RMatType>(ml,mr)/nrm;
//
//    return out;
//}
//
//MPSBlockMat<IKey,Complex>
//TransverseIsingModel::RandCMPSMat(size_t ml0, size_t mr0, size_t Nsecs) const//, size_t Nsites=1)
//{
//    size_t ml=ml0/2,mr=mr0/2,Nsites=1;
//    srand(time(NULL));
//
//    MPSBlockMat<IKey,Complex> out(this->localdim_,Nsites);
//    IKey K1({0}),K2({1});
//
//    out[0][IMPSKey(K1,this->Plus(K1,this->keys_[0]))] = randn<CMatType>(ml,mr);
//    out[0][IMPSKey(K2,this->Plus(K2,this->keys_[0]))] = randn<CMatType>(ml,mr);
//
//    out[1][IMPSKey(K1,this->Plus(K1,this->keys_[1]))] = randn<CMatType>(ml,mr);
//    out[1][IMPSKey(K2,this->Plus(K2,this->keys_[1]))] = randn<CMatType>(ml,mr);
//
//    return out;
//}


/**< FERMI HUBBARD MODEL TYPE *******************************************************************/
void FHubModel::InitSyms()
{
    /// FHUB has two U(1) symmetries: particle numbers for down (1st) and up (2nd)
    keys_.reserve(4);
    keys_.push_back({-1,-1}); /// |0>
    keys_.push_back({1,-1}); /// |->
    keys_.push_back({-1,1}); /// |+>
    keys_.push_back({1,1}); /// |+->
}

IKey FHubModel::Plus(const IKey& lhs, const IKey& rhs) const { return (lhs.Plus(rhs));}
IKey FHubModel::Minus(const IKey& lhs, const IKey& rhs) const { return (lhs.Minus(rhs));}

//MPSBlockMat<IKey,Real>
//FHubModel::RandRMPSMat(size_t ml0, size_t mr0, size_t Nsecs) const//, size_t Nsites=1)
//{
//    /**< Nsecs is the number of sectors for EACH quantum number! */
//    size_t Nsites=1;
//    size_t Nsecs2 = (Nsecs+1)*(Nsecs+1);
//    srand(time(NULL));
//
//    ivec secs;
//    if ( Nsecs%2 == 0) secs=linspace<ivec>(-Nsecs/2,Nsecs/2,Nsecs+1);
//    else secs = linspace<ivec>(-(int)(Nsecs/2+1),Nsecs/2,Nsecs+1);
//
//    /// uniform distribution of matrix dimensions into symmetry sectors
//    imat ml = (ml0/Nsecs2)*ones<imat>(Nsecs+1,Nsecs+1);
//    imat mr = (mr0/Nsecs2)*ones<imat>(Nsecs+1,Nsecs+1);
//
////    secs.print("secs");
////
////    cout<<ml0<<", "<<sum(sum(ml))<<endl;
////    cout<<mr0<<", "<<sum(sum(mr))<<endl;
////    double nrm=1;
//    double nrm=sqrt(2*sqrt(mean(vectorise(ml))*mean(vectorise(mr))));
////
//    MPSBlockMat<IKey,Real> out(this->localdim_,Nsites);
//    IKey K;
//    IMPSKey MK;
//    size_t s;
//
//    for (uint i=0; i<Nsecs; ++i)
//    {
//        for (uint j=0; j<Nsecs; ++j)
//        {
//            /// |0>
//            s = 0;
//            K = IKey({secs(i+1),secs(j+1)});
//            MK = MPSKeyType<IKey>(K,this->Plus(K,this->keys_[s]));
//            out[s][MK] = randn(ml(i+1,j+1),mr(i,j))/nrm;
//
//            /// |1> = |up>
//            s = 1;
//            K = IKey({secs(i),secs(j+1)});
//            MK = MPSKeyType<IKey>(K,this->Plus(K,this->keys_[s]));
//            out[s][MK] = randn(ml(i,j+1),mr(i+1,j))/nrm;
//
//
//            /// |2> = |down>
//            s = 2;
//            K = IKey({secs(i+1),secs(j)});
//            MK = MPSKeyType<IKey>(K,this->Plus(K,this->keys_[s]));
//            out[s][MK] = randn(ml(i+1,j),mr(i,j+1))/nrm;
//
//            /// |3> = |down,up>
//            s = 3;
//            K = IKey({secs(i),secs(j)});
//            MK = MPSKeyType<IKey>(K,this->Plus(K,this->keys_[s]));
//            out[s][MK] = randn(ml(i,j),mr(i+1,j+1))/nrm;
//        }
//    }
//    return out;
//}

/**< TEST 2*Z2 MODEL TYPE **********************************/
void Z2TestModel::InitSyms()
{
    keys_.reserve(4);
    keys_.push_back({0,0});
    keys_.push_back({1,0});
    keys_.push_back({0,1});
    keys_.push_back({1,1});
}


/**< For Z2, group addition and subtraction are the same! */
IKey Z2TestModel::Plus(const IKey& lhs, const IKey& rhs) const { return (lhs.Plus(rhs)).Modulus(2);}
IKey Z2TestModel::Minus(const IKey& lhs, const IKey& rhs) const { return Plus(lhs,rhs);} /// for
//IKey Z2TestModel::Plus(const IKey& lhs, const IKey& rhs) const
//
//IKey Z2TestModel::Minus(const IKey& lhs, const IKey& rhs) const

#endif // _USE_SYMMETRIES_


