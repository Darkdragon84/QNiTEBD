#ifndef MODELS_H_
#define MODELS_H_

//#warning "models included"

#include <iostream>
#include "arma_typedefs.h"
//#include "helpers.hpp"
#include "OperatorTypes.hpp"

#ifdef _USE_SYMMETRIES_
#include <map>
//#include "KeyTypes.hpp"
//#include "MPSBlockMat.hpp"
#endif // _USE_SYMMETRIES_

//#include "fulldiag.h"
template<typename VT>
class VecKeyType;

typedef VecKeyType<int> IKey;

/// GLOBAL HELPER FUNCTIONS FOR MODELS -----------------------------------------------------------/
const std::vector<size_t> num2ditvec(size_t x, size_t d, size_t N);
size_t ditvec2num(const std::vector<size_t>& vec,size_t d);

/**< BASE CLASS MODEL ******************************************************/
class ModelType
{
public:
    ModelType():localdim_(0),Nsym_(0) {};
    ModelType(size_t localdim, int Nsym):localdim_(localdim),Nsym_(Nsym) {};
    virtual ~ModelType() {};

    virtual void ShowParams() const = 0; /// define as pure virtual function
//    virtual const IKey Add(const IKey& lhs,const IKey& rhs) const = 0;

    inline const SparseOperator<double> Get2SiteHam() const { return ham_;};
    inline size_t GetLocalDim() const{return localdim_;};

    const SparseOperator<double> ImagTimeEvoOp(Real dt) const;

#ifdef _USE_SYMMETRIES_
    typedef typename std::vector<std::vector<int> > iivec;
//    inline void ShowSyms() const {size_t i=0;for (const auto& kit : keys_) cout<<i++<<" = "<<kit<<endl;}
    inline uint GetNsym() const {return Nsym_;};
//    virtual MPSBlockMat<IKey,Real> RandRMPSMat(size_t ml0, size_t mr0, size_t Nsecs) const = 0;
//    virtual MPSBlockMat<IKey,Complex> RandCMPSMat(size_t ml0, size_t mr0, size_t Nsecs) const = 0;

    /// GROUP OPERATIONS
    virtual IKey Plus(const IKey& lhs, const IKey& rhs) const = 0; /// define as pure virtual function
    virtual IKey Minus(const IKey& lhs, const IKey& rhs) const = 0; /// define as pure virtual function

    inline const std::vector<int>& GetRKey(size_t i) const {assert(i<localdim_); return keys_[i];};
    inline std::vector<int> GetKey(size_t i) const {assert(i<localdim_); return keys_[i];};
//    virtual const std::vector<IKey> IToKeyVec(const size_t Nsites) const = 0; /// define as pure virtual function
#endif // _USE_SYMMETRIES_
//    const std::vector<int>
protected:
    const size_t localdim_;
    const int Nsym_;
    SparseOperator<double> ham_,id_;
    virtual void MakeOps() = 0;

#ifdef _USE_SYMMETRIES_
    virtual void InitSyms() = 0;
//    std::vector<IKey> keys_;
    iivec keys_;
#endif // _USE_SYMMETRIES_
};

/**< XYZ MODEL BASE TYPE *******************************************************/
class XYZModelType : public ModelType
{
public:
    XYZModelType(double Jx,double Jy,double Jz,double hz,size_t dim,size_t Nsyms=1);
    virtual ~XYZModelType(){};
//    XYZModelType(double Jx,double Jy,double Jz,double hx,double hz,size_t dim);
    virtual void ShowParams () const;
    inline SparseOperator<Real> GetSx() {return sx_;};
    inline SparseOperator<Real> GetSyi() {return syi_;};
    inline SparseOperator<Real> GetSz() {return sz_;};

    #ifdef _USE_SYMMETRIES_
    virtual IKey Plus(const IKey& lhs, const IKey& rhs) const = 0; /// define as pure virtual function
    virtual IKey Minus(const IKey& lhs, const IKey& rhs) const = 0; /// define as pure virtual function
//    virtual const std::vector<IKey> IToKeyVec(const size_t Nsites) const = 0; /// define as pure virtual function
//    virtual MPSBlockMat<IKey,Real> RandRMPSMat(size_t ml0, size_t mr0, size_t Nsecs) const = 0; /// define as pure virtual function
//    virtual MPSBlockMat<IKey,Complex> RandCMPSMat(size_t ml0, size_t mr0, size_t Nsecs) const = 0; /// define as pure virtual function
    #endif // _USE_SYMMETRIES_
//    virtual const std::vector<IKey>& IToKeyVec(const size_t Nsites) const {};
//    const SparseOperator<double> Get2SiteHam() const;
protected:
    double Jx_,Jy_,Jz_,hz_;
    SparseOperator<double> sx_,syi_,sz_;
    std::map<uint,std::string> spin_;
    void MakeOps();
//#ifdef _USE_SYMMETRIES_
//    virtual void InitSyms() = 0;
//#endif // _USE_SYMMETRIES_
};


/**< XYZ MODEL TYPE *******************************************************/
class XYZModel : public XYZModelType
{
public:
    XYZModel(double Jx, double Jy, double Jz, double hz,size_t dim):XYZModelType(Jx,Jy,Jz,hz,dim,1)
    {
#ifdef _USE_SYMMETRIES_
        InitSyms();
        #endif // _USE_SYMMETRIES_
    };
    void ShowParams() const
    {
        cout<<spin_.find(localdim_)->second<<"XYZ type model:"<<endl;
        this->XYZModelType::ShowParams();
    }
#ifdef _USE_SYMMETRIES_
    IKey Plus(const IKey& lhs, const IKey& rhs) const;
    IKey Minus(const IKey& lhs, const IKey& rhs) const;
//    const std::vector<IKey> IToKeyVec(const size_t Nsites) const;
//    MPSBlockMat<IKey,Real> RandRMPSMat(size_t ml0, size_t mr0, size_t Nsecs) const;// {return MPSBlockMat<IKey,Real>(this->localdim_,1);};
//    MPSBlockMat<IKey,Complex> RandCMPSMat(size_t ml0, size_t mr0, size_t Nsecs) const;// {return MPSBlockMat<IKey,Complex>(this->localdim_,1);};
#endif // _USE_SYMMETRIES_
//    const IKey Add(const IKey& lhs,const IKey& rhs) const {};

protected:
#ifdef _USE_SYMMETRIES_
    void InitSyms();
#endif // _USE_SYMMETRIES_
};

//class XYZModel : public XYZModelType
//{
//public:
//    XYZModel(double Jx,double Jy,double Jz,double hz,size_t dim):XYZModelType(Jx,Jy,Jz,hz,dim,1)
//    {
//#ifdef _USE_SYMMETRIES_
//        InitSyms();
//        #endif // _USE_SYMMETRIES_
//    }
////    XYZModelType(double Jx,double Jy,double Jz,double hx,double hz,size_t dim);
//    void ShowParams () const;
//
//#ifdef _USE_SYMMETRIES_
//    const std::vector<IKey> IToKeyVec(const size_t Nsites) const;
//#endif // _USE_SYMMETRIES_
////    const IKey Add(const IKey& lhs,const IKey& rhs) const {};
////    const SparseOperator<double> Get2SiteHam() const;
//protected:
//#ifdef _USE_SYMMETRIES_
//    void InitSyms();
//#endif // _USE_SYMMETRIES_
//};

/**< XXZ MODEL TYPE *******************************************************/
class XXZModel : public XYZModelType
{
public:
    XXZModel(double J,double Delta,double hz,size_t dim):XYZModelType(J,J,Delta,hz,dim,1),J_(J),Delta_(Delta)
    {
#ifdef _USE_SYMMETRIES_
        InitSyms();
        #endif // _USE_SYMMETRIES_
    };
    void ShowParams() const
    {
        cout<<spin_.find(localdim_)->second<<"XXZ type model:"<<endl;
        this->XYZModelType::ShowParams();
    }
#ifdef _USE_SYMMETRIES_
    IKey Plus(const IKey& lhs, const IKey& rhs) const;
    IKey Minus(const IKey& lhs, const IKey& rhs) const;
//    const std::vector<IKey> IToKeyVec(const size_t Nsites) const;
//    MPSBlockMat<IKey,Real> RandRMPSMat(size_t ml0, size_t mr0, size_t Nsecs) const;// {return MPSBlockMat<IKey,Real>(this->localdim_,1);};
//    MPSBlockMat<IKey,Complex> RandCMPSMat(size_t ml0, size_t mr0, size_t Nsecs) const;// {return MPSBlockMat<IKey,Complex>(this->localdim_,1);};
#endif // _USE_SYMMETRIES_
//    const IKey Add(const IKey& lhs,const IKey& rhs) const {};

protected:
    double J_,Delta_;
#ifdef _USE_SYMMETRIES_
    void InitSyms();
#endif // _USE_SYMMETRIES_
};

/**< TRANSVERSE ISING MODEL TYPE *******************************************************/
class TransverseIsingModel : public XYZModelType
{
public:
    TransverseIsingModel(double J, double h):XYZModelType(J,0.,0.,h,2,1),J_(J),h_(h)
    {
#ifdef _USE_SYMMETRIES_
        InitSyms();
        #endif // _USE_SYMMETRIES_
    };
    void ShowParams() const
    {
        cout<<spin_.find(localdim_)->second<<"Transverse Ising type model:"<<endl;
        this->XYZModelType::ShowParams();
    }
#ifdef _USE_SYMMETRIES_
    IKey Plus(const IKey& lhs, const IKey& rhs) const;
    IKey Minus(const IKey& lhs, const IKey& rhs) const;
//    const std::vector<IKey> IToKeyVec(const size_t Nsites) const;
//    MPSBlockMat<IKey,Real> RandRMPSMat(size_t ml0, size_t mr0, size_t Nsecs=2) const;
//    MPSBlockMat<IKey,Complex> RandCMPSMat(size_t ml0, size_t mr0, size_t Nsecs=2) const;
#endif // _USE_SYMMETRIES_
//    const IKey Add(const IKey& lhs,const IKey& rhs) const{};

protected:
    double J_,h_;
#ifdef _USE_SYMMETRIES_
    void InitSyms();
#endif // _USE_SYMMETRIES_
};

/**< FERMI HUBBARD MODEL TYPE ********************************************/
class FHubModel : public ModelType
{
public:
    FHubModel(double t,double U,double mu,double hz);
    void ShowParams () const;
#ifdef _USE_SYMMETRIES_
    IKey Plus(const IKey& lhs, const IKey& rhs) const;
    IKey Minus(const IKey& lhs, const IKey& rhs) const;
//    const std::vector<IKey> IToKeyVec(const size_t Nsites) const;
//    MPSBlockMat<IKey,Real> RandRMPSMat(size_t ml0, size_t mr0, size_t Nsecs) const;// {return MPSBlockMat<IKey,Real>(this->localdim_,1);};
//    MPSBlockMat<IKey,Complex> RandCMPSMat(size_t ml0, size_t mr0, size_t Nsecs) const {return MPSBlockMat<IKey,Complex>(this->localdim_,1);};
#endif // _USE_SYMMETRIES_
//    const IKey Add(const IKey& lhs,const IKey& rhs) const{};

protected:
    double t_,U_,mu_,hz_;
    SparseOperator<double> cup_,cdo_,nup_,ndo_;
    void MakeOps();
#ifdef _USE_SYMMETRIES_
    void InitSyms();
#endif // _USE_SYMMETRIES_
};

/**< EXTENDED FERMI HUBBARD MODEL TYPE **********************************/
class ExtFHubModel : public FHubModel
{
public:
    ExtFHubModel(double t,double U,double mu,double hz,double V, double a);
    void ShowParams () const;
//    const SparseOperator<double> Get2SiteHam() const {};
protected:
    double V_,a_;
#ifdef _USE_SYMMETRIES_
    void InitSyms()
    {
        cerr<<"ExtHub::InitSyms() not implemented"<<endl;
    };
#endif // _USE_SYMMETRIES_
//    void MakeOps() {};
};

/**< TEST 2*Z2 MODEL TYPE **********************************/
class Z2TestModel : public ModelType
{
public:
    Z2TestModel();
    void ShowParams () const {};
#ifdef _USE_SYMMETRIES_
    IKey Plus(const IKey& lhs, const IKey& rhs) const;
    IKey Minus(const IKey& lhs, const IKey& rhs) const;
#endif // _USE_SYMMETRIES_
//    const SparseOperator<double> Get2SiteHam() const {};
protected:
    void MakeOps() {};
#ifdef _USE_SYMMETRIES_
    void InitSyms();
#endif // _USE_SYMMETRIES_
};

//#ifdef _USE_SYMMETRIES_
/////**< TEMPLATE CLASS ITOKEY TO GAIN INTERFACE BETWEEN KEYS AND MODELS ********/
//template<size_t Nsites,typename KT>
//class IToKey
//{
//public:
//    IToKey() {};
//    IToKey(const ModelType* model):model_(model),localdim_(model->GetLocalDim()),keys_(model->IToKeyVec(Nsites))
//    {
//        numel_=keys_.size();
//    };
//    ~IToKey() {};
//
//    inline KT Plus(const KT& in, size_t i) const {assert(i<numel_); KT out = model_->Plus(in,keys_[i]); return out;};
//    inline KT Minus(const KT& in, size_t i) const {assert(i<numel_); KT out = model_->Minus(in,keys_[i]); return out;};
//
//    inline size_t GetLocalDim() const {return localdim_;};
//    inline size_t GetNSites() const {return Nsites;};
//    inline size_t GetNumel() const {return numel_;};
//    inline const KT& operator[](size_t i) const /// return const reference
//    {
//        assert(i<numel_);
//        return keys_[i];
//    }
////    inline KT& operator[](size_t i) /// return non-const reference
////    {
////        assert(i<numel_);
////        return keys_[i];
////    }
//protected:
//    const ModelType* model_;
//    const size_t localdim_;
//    size_t numel_;
//    std::vector<KT> keys_;
//};


//#endif // _USE_SYMMETRIES_


#endif // MODELS_H_
