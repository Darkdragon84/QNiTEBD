#ifndef MPS_BLOCK_MAT_H
#define MPS_BLOCK_MAT_H

#include <vector>
#include <map>
#include <cmath>
#include <cassert>
#include <fstream>

#include "arma_typedefs.h"
#include "BlockObj.hpp"

using std::map;
using std::vector;
using std::cout;
using std::endl;



/** MPS MATRIX CLASS (WITH QUANTUM NUMBERS) *****************************************************************************************************
defined as a vector of maps. The vector indices correspond to the MPS physical indices and for each of them there is a matrix in blockform containing
values of VT type. These blockmatrices are stored as a map, where KT is the Key after which the separate blocks are ordered within the map
(see class BlockMat in BlockObj.hpp). In more detail, the map Key is actually the left (set of) quantum number(s), which is connected to the right (set of)
quantum number(s) via the physical index according to the underlying abelian group operation, i.e. [right = left + index] or [left = right - index].
a pair of KT, where the left and right element are connected via the physical index according to the underlying abelian group operation, i.e.
[right = left + index] or [left = right - index]. The blocks are therefore stored according to their left quantum numbers.

KT can be any type viable as a key for a std::map, i.e. it needs to have a copy (and if possible, move) constructor and there must be an operation < that compares
two elements of KT type. Additionally, for the MPS container, the + and - operators need to be defined. They should correspond to the abelian group action.
/// DECLARATION --------------------------------------------------------------------------------------------------------------------------------*/


template<typename KT,typename VT>
class MPSBlockMat : public std::vector<BlockMat<KT,VT> >
{
public:
    typedef typename BlockMat<KT,VT>::iterator miter;
    typedef typename BlockMat<KT,VT>::const_iterator mciter;

//    MPSBlockMat() = delete;
//    MPSBlockMat(){};
    MPSBlockMat(uint d,uint N):std::vector<BlockMat<KT,VT> >(pow(d,N)),d_(d),NSites_(N) {};

    inline uint GetLocalDim() const {return d_;};
    inline uint GetNSites() const {return NSites_;};

    /// helpers
    void ShowDims(const std::string& name="") const;
    map<KT,uint> GetAllSizes() const;
    inline uint GetTotalSizeLeft() const;
    inline uint GetTotalSizeRight() const;

    /// DiskIO
    bool Save(std::string name);
    bool Save(std::ofstream& file);
    bool Load(std::string name);
    bool Load(std::ifstream& file);

//    void MultLamLeft(const BlockLambda<KT>& lam);
//    void DivLamLeft(const BlockLambda<KT>& lam, double thresh=1e-15);
//    void MultLamRight(const BlockLambda<KT>& lam);
//    void DivLamRight(const BlockLambda<KT>& lam, double thresh=1e-15);

protected:
    uint d_,NSites_;
};

/**< std ostream output for e.g. screen output with cout */
template<typename KT, typename VT>
ostream&
operator<<(ostream& os, const MPSBlockMat<KT,VT>& MPS)
{
    for (uint i=0;i<MPS.size();++i)
    {
        os<<i<<":"<<endl;
        os<<MPS[i]<<endl;
    }
    return os;
}

/**< Show present quantum number sectors and corresponding matrix sizes for each physical index */
template<typename KT, typename VT>
void
MPSBlockMat<KT,VT>::ShowDims(const std::string& name) const
{
    if(name!="")cout<<name<<endl;
    for (uint i=0;i<this->size();++i)
    {
        cout<<i<<":"<<endl;
        this->at(i).ShowDims();
    }
}

template<typename KT,typename VT>
map<KT,uint>
MPSBlockMat<KT,VT>::GetAllSizes() const
{
    map<KT,uint> dimmap;
    for (const auto& Ait : *this)
    {
        for (const auto& matit : Ait)
        {
/// TODO (valentin#1#2015-05-12): consider using lower_bound to get rid of emplace()
            auto res1 = dimmap.emplace(Qin(matit),QMat(matit).n_rows);
            if ( !res1.second ) assert(res1.first->second == QMat(matit).n_rows && "GetAllSizes(): encountered same symmetry sector with different sizes") ;
            auto res2 = dimmap.emplace(Qout(matit),QMat(matit).n_cols);
            if ( !res2.second ) assert(res2.first->second == QMat(matit).n_cols && "GetAllSizes(): encountered same symmetry sector with different sizes") ;
        }
    }
    return dimmap;
}

template<typename KT,typename VT>
inline
uint
MPSBlockMat<KT,VT>::GetTotalSizeLeft() const
{
    uint dim=0;
    std::map<KT,uint> sizes;
    for (const auto& vit : *this) for(const auto& mit : vit) sizes.emplace(Qin(mit),QMat(mit).n_rows);
    for (const auto& sit : sizes) dim+= sit.second;
    return dim;
}

template<typename KT,typename VT>
inline
uint
MPSBlockMat<KT,VT>::GetTotalSizeRight() const
{
    uint dim=0;
    std::map<KT,uint> sizes;
    for (const auto& vit : *this) for(const auto& mit : vit) sizes.emplace(Qout(mit),QMat(mit).n_cols);
    for (const auto& sit : sizes) dim+= sit.second;
    return dim;
}

/**< MPSBLOCKMAT BASIC UTILITIES ***********************************************************************************************/

/**< concatenate two MPS matrices */
template<typename KT,typename VT>
MPSBlockMat<KT,VT> operator*(const MPSBlockMat<KT,VT>& lhs,const MPSBlockMat<KT,VT>& rhs)
{
    assert(lhs.GetLocalDim()==rhs.GetLocalDim() && "MPSBlockMat<KT,VT> operator*(): lhs and rhs have different LocalDim");
    MPSBlockMat<KT,VT> out(lhs.GetLocalDim(),lhs.GetNSites() + rhs.GetNSites());
/// TODO (valentin#1#2015-04-20): switch to iterators
    for (uint i=0;i<lhs.size();++i)
    {
        for (uint j=0;j<rhs.size();++j) out[i*rhs.size() + j]=lhs[i]*rhs[j];
    }
    return out;
}

/**< multiplication by BlockDiagMat from the left */
template<typename KT,typename VT>
inline
MPSBlockMat<KT,VT>
operator*(const BlockDiagMat<KT,VT>& mat, const MPSBlockMat<KT,VT>& MPSin)
{
    MPSBlockMat<KT,VT> MPSout(MPSin.GetLocalDim(),MPSin.GetNSites());
/// TODO (valentin#1#2015-04-19): switch to iterators
    for (uint i=0;i<MPSin.size();++i) MPSout[i] = mat*MPSin[i];

//    typename MPSBlockMat<KT,VT>::const_iterator init;
//    typename MPSBlockMat<KT,VT>::iterator outit = MPSout.begin();
//    for (init=MPSin.begin(); init!=MPSin.end(); init++,outit++) *outit = mat*(*init);

    return MPSout;
}

/**< multiplication by BlockDiagMat from the left */
template<typename KT,typename VT>
inline
MPSBlockMat<KT,VT>
operator*(const MPSBlockMat<KT,VT>& MPSin, const BlockDiagMat<KT,VT>& mat)
{
    MPSBlockMat<KT,VT> MPSout(MPSin.GetLocalDim(),MPSin.GetNSites());
/// TODO (valentin#1#2015-04-19): switch to iterators
    for (uint i=0;i<MPSin.size();++i) MPSout[i] = MPSin[i]*mat;

//    typename MPSBlockMat<KT,VT>::const_iterator init;
//    typename MPSBlockMat<KT,VT>::iterator outit = MPSout.begin();
//    for (init=MPSin.begin(); init!=MPSin.end(); init++,outit++) *outit = mat*(*init);

    return MPSout;
}

/**< multiplication and division by Schmidt values */
//template<typename KT,typename VT>
//inline MPSBlockMat<KT,VT>& MultMPSLamRight(MPSBlockMat<KT,VT>& MPS, const BlockLambda<KT>& lam) {for (auto& it : MPS) MultBMatLamRight(it,lam); return MPS;}
//
//template<typename KT,typename VT>
//inline MPSBlockMat<KT,VT>& MultMPSLamLeft(const BlockLambda<KT>& lam, MPSBlockMat<KT,VT>& MPS) {for (auto& it : MPS) MultBMatLamLeft(lam,it); return MPS;}

/// modifying
template<typename KT,typename VT>
inline MPSBlockMat<KT,VT>& operator<(MPSBlockMat<KT,VT>& MPS, const BlockLambda<KT>& lam) {for (auto& it : MPS) MultBlockMatLamRight(it,lam); return MPS;}

template<typename KT,typename VT>
inline MPSBlockMat<KT,VT>& operator>(const BlockLambda<KT>& lam, MPSBlockMat<KT,VT>& MPS) {for (auto& it : MPS) MultBlockMatLamLeft(lam,it); return MPS;}

template<typename KT,typename VT>
inline MPSBlockMat<KT,VT>& operator>(MPSBlockMat<KT,VT>& MPS, const BlockLambda<KT>& lam) {for (auto& it : MPS) DivBlockMatLamRight(it,lam); return MPS;}

template<typename KT,typename VT>
inline MPSBlockMat<KT,VT>& operator<(const BlockLambda<KT>& lam, MPSBlockMat<KT,VT>& MPS) {for (auto& it : MPS) DivBlockMatLamLeft(lam,it); return MPS;}

/// non-modifying
template<typename KT,typename VT>
inline MPSBlockMat<KT,VT> operator<<(const MPSBlockMat<KT,VT>& MPSin, const BlockLambda<KT>& lam) {MPSBlockMat<KT,VT> MPSout(MPSin); for (auto& it : MPSout) MultBlockMatLamRight(it,lam); return MPSout;}

template<typename KT,typename VT>
inline MPSBlockMat<KT,VT> operator>>(const BlockLambda<KT>& lam, const MPSBlockMat<KT,VT>& MPSin) {MPSBlockMat<KT,VT> MPSout(MPSin); for (auto& it : MPSout) MultBlockMatLamLeft(lam,it); return MPSout;}

template<typename KT,typename VT>
inline MPSBlockMat<KT,VT> operator>>(const MPSBlockMat<KT,VT>& MPSin, const BlockLambda<KT>& lam) {MPSBlockMat<KT,VT> MPSout(MPSin); for (auto& it : MPSout) DivBlockMatLamRight(it,lam); return MPSout;}

template<typename KT,typename VT>
inline MPSBlockMat<KT,VT> operator<<(const BlockLambda<KT>& lam, const MPSBlockMat<KT,VT>& MPSin) {MPSBlockMat<KT,VT> MPSout(MPSin); for (auto& it : MPSout) DivBlockMatLamLeft(lam,it); return MPSout;}

/**< DISK IO ******************************************************************************************/
template<typename KT,typename VT>
bool
MPSBlockMat<KT,VT>::Save(std::string name)
{
//    std::string fullname = name + ".bin";
    std::ofstream file(name.c_str(), std::fstream::binary);
    bool save_okay = Save(file);
    file.close();

    return save_okay;
}

template<typename KT,typename VT>
bool
MPSBlockMat<KT,VT>::Save(std::ofstream& file)
{
    bool save_okay = file.is_open();
    uint Nsec=0;

    file.write(reinterpret_cast<char*>(&d_), sizeof(uint));
    file.write(reinterpret_cast<char*>(&NSites_), sizeof(uint));

    for (uint s=0;s<this->size();++s)
    {
        Nsec=this->at(s).size();
        file.write(reinterpret_cast<char*>(&Nsec), sizeof(uint));
        for (const auto& it : this->at(s))
        {
            file << Qin(it);
            file << Qout(it);
            save_okay = save_okay && diskio::save_arma_binary(QMat(it),file);
        }
    }
    return save_okay;
}

template<typename KT,typename VT>
bool
MPSBlockMat<KT,VT>::Load(std::string name)
{
//    std::string fullname = name + ".bin";
    this->clear();
    std::ifstream file(name.c_str(), std::fstream::binary);
    bool save_okay = file.is_open();
    uint Nsec=0;

    file.read(reinterpret_cast<char*>(&d_), sizeof(uint));
    file.read(reinterpret_cast<char*>(&NSites_), sizeof(uint));
    this->resize(pow(d_,NSites_));

    Mat<VT> tmp;
    std::string err;

    KT Kin,Kout;
    for (auto& Ait : *this)
    {
        file.read(reinterpret_cast<char*>(&Nsec), sizeof(uint));
        for (uint i=0;i<Nsec;++i)
        {
            file>>Kin;
            file>>Kout;
            diskio::load_arma_binary(tmp,file,err);

            Ait.emplace_hint(Ait.end(),Kin,QMatPair<KT,VT>(Kout,tmp));
        }
    }

    file.close();

    return save_okay;
}


#endif // MPS_BLOCK_MAT_H
