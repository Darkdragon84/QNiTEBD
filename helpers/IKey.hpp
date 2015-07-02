#ifndef IKEY_H_
#define IKEY_H_

#include <vector>
#include <cassert>
#include <algorithm>
#include <initializer_list>

//class IKey : public std::vector<int>
//{
//public:
//
//};
typedef std::vector<int> IKey;

//IKey operator+(const IKey& lhs, const IKey& rhs)
//{
//    IKey out(lhs.size());
//    std::transform(lhs.begin(),lhs.end(),rhs.begin(),out.begin(),std::plus<int>());
//    return out;
//}
//
//IKey operator-(const IKey& lhs, const IKey& rhs)
//{
//    IKey out(lhs.size());
//    std::transform(lhs.begin(),lhs.end(),rhs.begin(),out.begin(),std::minus<int>());
//    return out;
//}

using std::ostream;
using std::ofstream;
using std::ifstream;

IKey operator+(const IKey& lhs, const IKey& rhs)
{
    uint siz=lhs.size();
    IKey out(siz);
    if (siz==3) {for (uint i=1;i<3;++i)out[i]=lhs[i]+rhs[i];out[0]=(lhs[0]+rhs[0])%2;}
    else for (uint i=0;i<siz;++i)out[i]=lhs[i]+rhs[i];
    return out;
}

IKey operator-(const IKey& lhs, const IKey& rhs)
{
    uint siz=lhs.size();
    IKey out(siz);
    if (siz==3) {for (uint i=1;i<3;++i)out[i]=lhs[i]-rhs[i];out[0]=(lhs[0]+rhs[0])%2;}
    else for (uint i=0;i<siz;++i)out[i]=lhs[i]-rhs[i];
    return out;
}


ostream& operator<<(ostream& os,const IKey& K)
{
    for (unsigned int i=0;i<K.size();++i) os<<(i==0 ? "(":",")<<K[i];
    os<<")";
    return os;
}

ofstream& operator<<(ofstream& file, const IKey& K)
{
    unsigned int Nsym=K.size();
//    file<<Nsym;
//    for (unsigned int i=0;i<Nsym;++i) file<<K[i];
//    int QN=0;
    file.write(reinterpret_cast<char*>(&Nsym),sizeof(unsigned int));
    for (unsigned int i=0;i<Nsym;++i) file.write(reinterpret_cast<const char*>(&K[i]),sizeof(int));
//    for (unsigned int i=0;i<K.size();++i) file<<(i>0 ? " ":"")<<K[i];
    return file;
}

ifstream& operator>>(ifstream& file, IKey& K)
{
    unsigned int Nsym=0;
    file.read(reinterpret_cast<char*>(&Nsym),sizeof(unsigned int));
    K.resize(Nsym);

    for (unsigned int i=0;i<Nsym;++i) file.read(reinterpret_cast<char*>(&K[i]),sizeof(int));
    return file;
}

template<uint N_>
class ItoKey : public std::vector<IKey>
{
public:
//    ItoKey(){};
    ItoKey(uint d): d_(d) {};
    ItoKey(uint d, const std::initializer_list<IKey>& lst):vector<IKey>(lst),d_(d) {assert(lst.size()==pow(d_,N_));};

    inline uint GetLocalDim() const {return d_;};
    inline uint GetNSites() const {return N_;};
protected:
    uint d_;
};

#endif // IKEY_H_
