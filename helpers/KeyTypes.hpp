#ifndef KEY_TYPES_H
#define KEY_TYPES_H

//#warning "keytypes included"
#include <vector>
#include <assert.h>
#include <iostream>
#include <stdlib.h>
#include <initializer_list>
#include <functional>
#include <algorithm>
#include <utility>
#include <cmath>

#include "Models.h"

using namespace std;

/// *** Vector Key Type ************************************************************************************************************************//
template<typename VT>
class VecKeyType : public std::vector<VT>
{
public:
/// TODO (valentin#1#01/09/2014): Check if constructor for vector is really necessary or if it is enough to just reserve
//    VecKeyType() = default;
    VecKeyType(const ModelType& mod):std::vector<VT>(mod.GetNsym()),mod_(mod) {DOUT("KEY std");};

    VecKeyType(const VecKeyType& rhs):std::vector<VT>(rhs),mod_(rhs.GetRMod()) {DOUT("KEY copied");}; /// redefine copy ctor as move assignment definition breaks synthesized one
    VecKeyType(VecKeyType&& rhs):std::vector<VT>(std::move(rhs)),mod_(rhs.GetRMod()) {DOUT("KEY moved");}; /// move ctor as move assignment definition breaks synthesized one

/// TODO (valentin#1#01/09/2014): replace asserts by proper checks with error handling
    explicit VecKeyType(const ModelType& mod, const std::initializer_list<VT>& lst):std::vector<VT>(lst),mod_(mod)
    {
        assert(lst.size()==mod.GetNsym());
        DOUT("KEY created");
    };
    explicit VecKeyType(const ModelType& mod, std::initializer_list<VT>&& lst):std::vector<VT>(std::move(lst)),mod_(mod)
    {
        assert(this->size()==mod.GetNsym());
        DOUT("KEY created");
    };

    explicit VecKeyType(const ModelType& mod, const std::vector<VT>& vec):std::vector<VT>(vec),mod_(mod)
    {
        assert(this->size()==mod.GetNsym());
        DOUT("KEY created");
    };
    explicit VecKeyType(const ModelType& mod, std::vector<VT>&& vec):std::vector<VT>(std::move(vec)),mod_(mod)
    {
        assert(this->size()==mod.GetNsym());
        DOUT("KEY created");
    };


    VecKeyType& operator=(const VecKeyType& rhs);/// for now each KeyType object can only be associated with one and the same model object for its entire lifetime! is this necessary?
    /// define move assignment operator for refs. to rvalue (i.e. temporaries such as results from arithmetics)
    VecKeyType& operator=(VecKeyType&& rhs);/// for now each KeyType object can only be associated with one and the same model object for its entire lifetime! is this necessary?

    /// ModelType group operations as friends to access low level arithmetics
//    friend IKey ModelType::Plus(const IKey& lhs, const IKey& rhs) const;
//    friend IKey ModelType::Minus(const IKey& lhs, const IKey& rhs) const;

    inline const ModelType& GetRMod() const
    {
        return mod_;
    };


    /// low level arithmetics
    inline VecKeyType Plus(const VecKeyType& rhs) const;
    inline VecKeyType Minus(const VecKeyType& rhs) const;
    inline VecKeyType Modulus(const VT mod) const;
//    ~VecKeyType() {};
protected:

    const ModelType& mod_; /// This causes VeckKeyType not to have a default constructor!!
};

/// for assignments mod_ doesn't need to be assigned, as it is required to already be the same anyways (do the check though!)
template<typename VT>
VecKeyType<VT>& VecKeyType<VT>::operator=(const VecKeyType& rhs) /// std implementation
{
    assert(&mod_==&rhs.GetRMod());
    if (this != &rhs)
    {
        VecKeyType tmp(rhs);
        this->swap(tmp);
    }
    DOUT("KEY copy assigned");
    return *this;
}

template<typename VT>
VecKeyType<VT>& VecKeyType<VT>::operator=(VecKeyType&& rhs)
{
    assert(&mod_==&rhs.GetRMod());/// nothing to move here, they have to be the same anyways
    if (this != &rhs) std::vector<VT>::operator=(std::move(rhs)); /// use the std::vector move assignment operator to move the vector part of the key
    DOUT("KEY move assigned");
//    if (this != &rhs)this->swap(rhs);
    return *this;
}

/// *** Tensor Key Type as pair of arbitrary Key Types ********************************************************************************************//
/** TODO (valentin#1#09/19/2013): It is possible to remove the explicit keyword for the MPSKeyType constructor taking a single KeyType argument.
This way implicit type conversion to a dummy MPSKeyType containing the same key twice is possible. This would be useful to simplify finding single
keys in MPSBlockMats, but it is also dangerous to be accidentally used somewhere else. Leave it like this for now... */

//template<typename KT>
//class MPSKeyType : public std::pair<KT,KT>
//{
//public:
//    MPSKeyType() = delete; /**< no default constructor */
//    MPSKeyType(const ModelType& mod):std::pair<KT,KT>(KT(mod),KT(mod)),mod_(mod) {}; /**< extended std ctor */
//    MPSKeyType(const KT& lhs, const KT& rhs):std::pair<KT,KT>(lhs,rhs),mod_(lhs.GetRMod()) /**< copy ctor */
//    {
//        assert(&lhs.GetRMod()==&rhs.GetRMod());
//    };
//    MPSKeyType(KT&& lhs, KT&& rhs):std::pair<KT,KT>(std::move(lhs),std::move(rhs)),mod_(this->first.GetRMod())  /**< move ctor for both */
//    {
//        assert(&this->first.GetRMod()==&this->second.GetRMod());
//    };
//    MPSKeyType(const KT& lhs, KT&& rhs):std::pair<KT,KT>(lhs,std::move(rhs)),mod_(this->first.GetRMod())  /**< copy ctor for left, move ctor for right */
//    {
//        assert(&this->first.GetRMod()==&this->second.GetRMod());
//    };
//    MPSKeyType(KT&& lhs, const KT& rhs):std::pair<KT,KT>(std::move(lhs),rhs),mod_(this->first.GetRMod())  /**< move ctor for right, copy ctor for left */
//    {
//        assert(&this->first.GetRMod()==&this->second.GetRMod());
//    };
////    explicit MPSKeyType(const KT& lhs):std::pair<KT,KT>(lhs,lhs) {}; /// dummy constructor to mimic a single key. explicit keyword: see TODO item above
//protected:
//    const ModelType& mod_;
//};

/// typedefs
typedef VecKeyType<int> IKey;

/// NON-MODIFYING OPERATIONS ----------------------------------------------------------------------/

template<typename VT>
inline VecKeyType<VT> VecKeyType<VT>::Plus(const VecKeyType<VT>& rhs) const
{
    DOUT("plus start");
    assert(&mod_==&rhs.GetRMod());
/// TODO (valentin#1#01/08/2014): size check should be obsolete and covered by dimensionality of model
    assert(this->size()==rhs.size());
    VecKeyType<VT> out(mod_);
    std::transform(this->begin(),this->end(),rhs.begin(),out.begin(),std::plus<VT>());
    DOUT("plus end");
    return out;
}


template<typename VT>
inline VecKeyType<VT> VecKeyType<VT>::Minus(const VecKeyType<VT>& rhs) const
{
    assert(&mod_==&rhs.GetRMod());
/// TODO (valentin#1#01/08/2014): size check should be obsolete and covered by dimensionality of model
    assert(this->size()==rhs.size());
    VecKeyType<VT> out(mod_);
    std::transform(this->begin(),this->end(),rhs.begin(),out.begin(),std::minus<VT>());
    return out;
}

template<typename VT>
inline VecKeyType<VT> VecKeyType<VT>::Modulus(const VT modulo) const
{
    VecKeyType<VT> out(mod_);
//    std::transform(this->begin(),this->end(),out.begin(),UnaryModulus<VT>(modulo));
    std::transform(this->begin(),this->end(),out.begin(),[&](const VT& in){return in % modulo;}); /**< define "unary" modulus as a lambda function */
    return out;
}

/// GROUP OPERATION OPERATORS DEPENDENT ON MODELTYPE -------------------------------------------------------------------------------/
template<typename T>
inline VecKeyType<T> operator+(const VecKeyType<T>& lhs, const VecKeyType<T>& rhs)
{
    return lhs.GetRMod().Plus(lhs,rhs);
}

template<typename T>
inline VecKeyType<T> operator-(const VecKeyType<T>& lhs, const VecKeyType<T>& rhs)
{
    return lhs.GetRMod().Minus(lhs,rhs);
}

/// OUTPUT ---------------------------------------------------------------------------------------------------------------------------/
template<typename T>
ostream& operator<<(ostream& os, const VecKeyType<T>& vec) /**< to be able to pass KeyTypes to output streams (such as cout) */
{
    bool first=true;
    os<<"(";
    for (const auto& it:vec)
    {
        first ? first=false : os<<",";
        os<<it;
    }
    os<<")";
    return os;
}


///*** TEMPLATE CLASS ITOKEY TO GAIN INTERFACE BETWEEN KEYS AND MODELS **********************************************************************************//
template<size_t Nsites,typename KT> /// KT = KeyType
class IToKey : public std::vector<KT>
{
public:
//    IToKey() {};
    IToKey(const ModelType& model):model_(model),localdim_(model.GetLocalDim()) /**< call standard constructor of std::vector<KT> with 0 elements as KT has no default constructor */
    {
        CreateKeys();
    };
    ~IToKey() {};

    inline size_t GetLocalDim() const
    {
        return localdim_;
    };
    inline size_t GetNSites() const
    {
        return Nsites;
    };
protected:
    void CreateKeys();

    const ModelType& model_; /**< This causes IToKey not to have a default constructor */
    const size_t localdim_;
};

template<size_t Nsites,typename KT>
void IToKey<Nsites,KT>::CreateKeys()
{
    std::vector<KT> keys1;
    keys1.reserve(localdim_); /**< to avoid copying */


/**<
- an rvalue is being passed, so the move constructor for the Key will be called
- GetKey only passes a vector of ints, not a key, thus we have to construct a dummy key from this vector of ints
*/
    for (size_t i=0; i<localdim_; ++i) keys1.push_back(KT(model_,model_.GetKey(i)));


    if (Nsites==1) this->swap(keys1);
    else
    {
        size_t numel = pow(localdim_,Nsites);
        this->reserve(numel);/**< to avoid copying */

        std::vector<size_t> inds; /**< vector of local indices for each site */
        for (size_t i=0; i<numel; ++i)
        {
            inds = num2ditvec(i,localdim_,Nsites);
            KT tmp(model_,keys1[inds[0]]);

            for (size_t j=1; j<Nsites; ++j) tmp = tmp + keys1[inds[j]];

/// TODO (valentin#1#2015-03-03): Check if it is safe to move tmp into the back of IToKey, for now just copy (only done once anyways)
            this->push_back(std::move(tmp));
//            this->push_back(tmp);

        }
    }
}

/**< OUTPUT */
template<size_t Nsites,typename KT>
ostream& operator<<(ostream& os, const IToKey<Nsites,KT>& I2K)
{
    for (const auto& it:I2K)os<<it<<endl;
    return os;
}

/**< partial template specialization for IKeys (C++11) */
template<size_t Nsites>
using IToIKey = IToKey<Nsites,IKey>;

#endif // KEY_TYPES_H
