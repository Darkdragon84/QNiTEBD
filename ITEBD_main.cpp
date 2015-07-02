#include <iostream>

#include "../../helpers/arma_typedefs.h"
#include "../../helpers/parser.h"

#include "../../helpers/BlockObj.hpp"
#include "../../helpers/EigsBlockFunctions.hpp"
#include "../../helpers/helpers.hpp"
#include "../../helpers/MPSBlockUtilities.hpp"
#include "../../helpers/OperatorTypes.hpp"
#include "../../helpers/IKey.hpp"
#include "../../helpers/tictoc.hpp"

using namespace std;

void makei2k(ItoKey<1>& I2K, ItoKey<2>& I2K2, uint d);
map<IKey,uint> makedims(ItoKey<1>& I2K, ItoKey<2>& I2K2,uint d, uint m0, uint Ns);
std::pair<RSpOp,RSpOp> generateOps(uint d);

typedef MPSBlockMat<IKey,Real> IMPSBRM;
typedef BlockDiagMat<IKey,Real> IBDiagRM;
typedef BlockLambda<IKey> IBLam;

int main(int argc, char** argv)
{
    uint d=2,m0=2,Ns=5,mmax=200;
    Real lamthresh=1e-10,dt=1e-2,dtmin=1e-5,thresh=1e-10;
    std::vector<uint> mmaxv;
    parser pp(argc,argv);
    pp.GetValue(d,"d");
    pp.GetValue(m0,"m0");
//    pp.GetValue(mmax,"mmax");
    pp.GetValue(mmaxv,"mmaxv");
    pp.GetValue(Ns,"Ns");
    pp.GetValue(dtmin,"dt");
    pp.GetValue(lamthresh,"lamthresh");
    pp.GetValue(thresh,"thresh");

    cout.precision(10);

    ItoKey<1> I2K(d);
    ItoKey<2> I2K2(d);
    makei2k(I2K,I2K2,d);
    auto dimmap = makedims(I2K,I2K2,d,m0,Ns);

    auto ops = generateOps(d);
    RSpOp H = ops.first;
    RSpOp obsop = ops.second;
    RSpOp U = expmat(-dt*H);
    IMPSBRM A(d,1),B(d,1),psi(d,2);
    IBLam lamA,lamB;


    IMPSBRM A0 = RandUMPS<Real>(dimmap,I2K);

    A = OrthoUMPS(lamA,A0,c);
    B = A;
    lamB = lamA;
    psi = ((lamB>>A)<<lamA) * (B << lamB);
//    IMPSBRM B0 = RandUMPS<Real>(dimmap,I2K);
//    if(A0.Save("A.bin"))cout<<"A.bin saved"<<endl;
//    if(B0.Save("B.bin"))cout<<"B.bin saved"<<endl;
//    IMPSBRM AA = RandUMPS<Real>(dimmap,I2K2);
//    AA.ShowDims("AA");

//    for (const auto& it : dimmap) lamB.emplace_hint(lamB.end(),it.first,exp(randn(it.second,1)));
//    split(A0*B0,A,B,lamA,I2K);
//    lamA/=norm(lamA);
//    truncate(A,B,lamA,I2K,mmax,lamthresh);
//    psi = A*(lamA>>B);

//    cout<<"U"<<endl<<U;

    A.ShowDims("A");
    B.ShowDims("B");
    lamA.ShowDims("lamA");
    lamB.ShowDims("lamB");

//    Eold = trace(ApplyTMmixedLeft(ApplyOperator(psi,H,I2K2),psi));

//    for (uint ct=0;ct<100;++ct)
    Real Eold=0,Enew=0,E1=0,E2=0,dE=1;
    uint ct=0,measct=0,measint=100;
    Real nrm,obs1,obs2;
    IMPSBRM M1(d,1),M2(d,1);


    while (dt>dtmin)
    {
        for (uint n=0; n<mmaxv.size(); ++n)
        {
            mmax = mmaxv[n];
            cout<<"mmax="<<mmax<<endl;
            dE=1;Eold=0;Enew=0;
            while (abs(dE)>thresh*measint)
            {
                psi = ApplyOperator(psi,U,I2K2);
                split(psi,A,B,lamA,I2K);
                lamA/=norm(lamA);
                truncate(A,B,lamA,I2K,mmax,lamthresh);
                psi = (lamA>>(B>>lamB)) * (A<<lamA);

                if (measct==measint) E1 = trace(ApplyTMmixedLeft(ApplyOperator(psi,H,I2K2),psi))/trace(ApplyTMLeft(psi));

                psi = ApplyOperator(psi,U,I2K2);
                split(psi,B,A,lamB,I2K);
                lamB/=norm(lamB);
                truncate(B,A,lamB,I2K,mmax,lamthresh);
                psi = (lamB>>(A>>lamA)) * (B<<lamB);

                if (measct==measint)
                {
                    measct = 0;
                    nrm = trace(ApplyTMLeft(psi));
                    E2 = trace(ApplyTMmixedLeft(ApplyOperator(psi,H,I2K2),psi))/nrm;
                    M1 = B<<lamB;
                    M2 = lamB>>A;
                    obs1 = trace(ApplyTMmixedLeft(ApplyOperator(M1,obsop,I2K),M1))/trace(ApplyTMLeft(M1));
                    obs2 = trace(ApplyTMmixedLeft(ApplyOperator(M2,obsop,I2K),M2))/trace(ApplyTMLeft(M2));

                    Enew = 0.5*(E1 + E2);
                    dE = Enew - Eold;
//            if (abs(dE)<thresh)
//            {
//                dt*=0.1;
//                U = expmat(-dt*H);
//            }
                    Eold = Enew;
                    lamA.ShowDims("lamA");
                    lamB.ShowDims("lamB");
//            CheckOrtho(lamA<<B,lamB,c);
//            CheckOrtho(A>>lamA,lamB,c);
                    cout<<ct++<<": m="<<lamA.GetTotalSize()<<", dt="<<dt<<", E="<<Enew<<", dE="<<dE<<", nrm="<<nrm<<", obs=("<<obs1<<","<<obs2<<")"<<endl;
                }
                measct++;
            }
        }
        dt*=0.1;
        U = expmat(-dt*H);
        cout<<"dt="<<dt<<endl;
    }
    return 0;
}

std::pair<RSpOp,RSpOp> generateOps(uint d)
{
    RSpOp H(d,2),op1(d,1),op2(d,1),op3(d,1),obsop(d,1);
    RSpOp Id = RSpOp::SpId(d,1);
    switch (d)
    {
    case 2:
        op1(0,1)=0.5;
        op1(1,0)=0.5;

        op2(0,1)=-0.5;
        op2(1,0)=0.5;

        op3(0,0)=0.5;
        op3(1,1)=-0.5;

        obsop = op3;

        H = -kron(op1,op1) + kron(op2,op2) + kron(op3,op3);
        break;
    case 4:
        op1(0,2)=1.;
        op1(1,3)=1.;

        op2(0,1)=1.;
        op2(2,3)=-1.;

        op3(0,0)=1.;
        op3(1,1)=-1.;
        op3(2,2)=-1.;
        op3(3,3)=1.;

        obsop(1,1)=1.;
        obsop(2,2)=1.;
        obsop(3,3)=2.;

        H = -kron(op1.t()*op3,op1) + kron(op1*op3,op1.t()) - kron(op2.t()*op3,op2) + kron(op2*op3,op2.t());
        break;
    default:
        cerr<<"not implemented"<<endl;
    }
    return std::make_pair(H,obsop);
}

void makei2k(ItoKey<1>& I2K, ItoKey<2>& I2K2, uint d)
{
    switch (d)
    {
    case 2:
        I2K = ItoKey<1>(d, {{-1},{1}});
        break;
    case 3:
        I2K = ItoKey<1>(d, {{-1,-1},{-1,1},{1,-1}}); /// (nup=1/2, ndo=1/2)
        break;
    case 4:
        I2K = ItoKey<1>(d, {{-1,-1},{-1,1},{1,-1},{1,1}}); /// (nup=1/2, ndo=1/2)
//        I2K = ItoKey<1>(d,{{-2,-2},{-2,1},{1,-2},{1,1}}); /// (nup=2/3, ndo=2/3)
//        I2K = ItoKey<1>(d,{{-1,0},{0,-1},{0,1},{1,0}}); /// (n=1, m=0)
//        I2K = ItoKey<1>(d,{{-1,0},{1,-1},{1,1},{3,0}}); /// (n=1/2, m=0)
        break;
    default:
        cerr<<"d="<<d<<" not implemented"<<endl;
    }
    I2K2.clear();
    for (uint i=0; i<d; ++i) for (uint j=0; j<d; ++j) I2K2.push_back(I2K[i] + I2K[j]);
}

map<IKey,uint> makedims(ItoKey<1>& I2K, ItoKey<2>& I2K2,uint d, uint m0, uint Ns)
{
    map<IKey,uint> dimmap;
    vector<uint> sdims(Ns);
    for (uint i=0; i<ceil(Ns/2.); ++i)
    {
//        uint tmpd = pow(2,i+2);
        uint tmpd = (i+1)*m0;
        sdims[i] = tmpd;
        sdims[Ns-i-1] = tmpd;
    }

    switch (d)
    {
    case 2:
        for (uint i=0; i<Ns; ++i) dimmap.emplace_hint(dimmap.end(),IKey {int(i-Ns/2)},sdims[i]);
        break;
    case 3:
        for (uint i=0; i<Ns; ++i) for (uint j=0; j<Ns; ++j) dimmap.emplace_hint(dimmap.end(),IKey {int(i-Ns/2),int(j-Ns/2)},sdims[i]*sdims[j]);
        break;
    case 4:
        for (uint i=0; i<Ns; ++i) for (uint j=0; j<Ns; ++j) dimmap.emplace_hint(dimmap.end(),IKey {int(i-Ns/2),int(j-Ns/2)},sdims[i]*sdims[j]);
        break;
    default:
        cerr<<"d="<<d<<" not implemented"<<endl;
    }
    return dimmap;
}
