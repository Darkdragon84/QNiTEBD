
#include "eigs.h"

using namespace std;
//auto cc2fc = reinterpret_cast<double*>;

int eigs_cn(std::function<void (Complex*,Complex*)> MultOPx, int N, CVecType& vals, CMatType& vecs, int nev, std::string whch, double tol, int maxit, int ncv)
{
/// TODO (valentin#1#2015-06-18): Replace all dynamically allocated arrays with arma vectors
    vals.reset();
    vecs.reset();
    /// PARAMS FOR ZNAUPD ---------------------------------------------------------------------------------------------------------//
    int mode=1; /// standard EV problem A*x = lam*x
    maxit=std::min(maxit,N);

    int IDO=0;
    char BMAT='I';
    char WHICH[3];
    strcpy(WHICH,whch.c_str());
    int NEV=nev;
    double TOL=tol;
    int NCV = (ncv==0) ? std::min(std::max(20,2*NEV+1),N) : ncv; /// from MATLAB:eigs:872
    int LDV = N; /// even though the C++ equivalent of a COMPLEX*16 array of length N is a double array of length 2*N, LDV is still N, as fortran jumps by 2 indices through the C style double array.
    int LWORKL=3*NCV*NCV + 5*NCV;
    int IPARAM[11] = {1,0,maxit,1,0,0,mode,0,0,0,0};
    int IPNTR[14];
    int INFONAUP=1;

    CVecType RESID(N,fill::randn); /// on entry contains starting vector for Arnoldi, on exit contains the final residual vector
    CMatType V(N,NCV); /// will contain Krylov Basis from Arnoldi
    CVecType WORKD(3*N);
    CVecType WORKL(LWORKL);
    RVecType RWORK(NCV);
//    Complex * WORKD = new Complex[3*N];
//    Complex * WORKL = new Complex[LWORKL];
//    double * RWORK = new double[NCV];
    /// ZNAUPD -------------------------------------------------------------------------------------------------------------------------//
//    int dummy = 0;
    while (IDO!=99)
    {
        /// standard EV problem: A*x = lambda*x (here OP=A and B=I)
        /// generalized EV problem: A*x = lambda*M*x (here OP = inv[M]*A and B=M)
        znaupd_(&IDO,&BMAT,&N,WHICH,&NEV,&TOL,RESID.memptr(),&NCV,V.memptr(),&LDV,IPARAM,IPNTR,WORKD.memptr(),WORKL.memptr(),&LWORKL,RWORK.memptr(),&INFONAUP);

        switch (IDO)
        {
        case -1: /// compute Y = OP * X
            /// initialization phase, which is somehow never used...
            cerr<<"-1 not implemented"<<endl;
            break;
        case 1: /// compute Z = B * X and Y = OP * Z
            /// For standard EV Problem Y = Z as B = I (see above)
            MultOPx(&WORKD[IPNTR[0]-1],&WORKD[IPNTR[1]-1]);
            break;
        case 2: /// compute Y = M * X (only for generalized)
            cerr<<"2 not implemented"<<endl;
            break;
        case 3: /// calculate shifts, if specified to do this separately (if IPARAM(1) = 0)
            cerr<<"3 not implemented"<<endl;
            break;
        case 4: /// compute Z = OP * X, documentation doesn't even say, which memory to use, probably never used...
            cerr<<"4 not implemented"<<endl;
            break;
        case 99:/// ARPACK HAS CONVERGED
//            cout<<"convergence"<<endl;
            break;
        default:
            cerr<<"IDO has unknown value"<<endl;
        }
    }
    int nconv=0;

    if (IDO==99 && INFONAUP==0)
    {
        nconv=IPARAM[4];
//        cout<<"ARPACK: "<<nconv<<" eigenpairs have converged"<<endl;
//        cout<<IPARAM[8]<<" OPx"<<endl;
    }
    else
    {
        cerr<<"no convergence in ZNAUPD, on exit: "<<INFONAUP<<endl;

//        delete[] WORKD;
//        delete[] WORKL;
//        delete[] RWORK;
        return 0;
    }

    /// PARAMS FOR ZNEUPD -------------------------------------------------------------------------------------------------------------------------//
    int INFONEUP=0;
    int RVEC=1;
    char HOWMNY='A';
    IVecType SELECT(NCV,fill::zeros);
//    CVecType EV(NEV+1);
    vals.resize(NEV+1);
//    SELECT.zeros();
    vecs.resize(N,NEV+1);
//    CMatType Zmat(N,NEV+1);
//    Complex* D = EV.memptr();
//    Complex* Z = Zmat.memptr();
    int LDZ=N;
    Complex SIGMA(0,0);
    CVecType WORKEV(2*NCV);

//    DOUT("calling zneupd");
    zneupd_(&RVEC,&HOWMNY,SELECT.memptr(),vals.memptr(),vecs.memptr(),&LDZ,&SIGMA,WORKEV.memptr(),&BMAT,&N,WHICH,&NEV,&TOL,RESID.memptr(),
            &NCV,V.memptr(),&LDV,IPARAM,IPNTR,WORKD.memptr(),WORKL.memptr(),&LWORKL,RWORK.memptr(),&INFONEUP);
//    DOUT("done");
    if (INFONEUP==0)
    {
        vals.resize(nconv);
        vecs.resize(N,nconv);
//        cout<<"ZNEUPD converged:"<<endl;
//        cout<<vals<<endl;
    }
    else
    {
        cerr<<"no convergence in ZNEUPD, on exit "<<INFONEUP<<endl;

//        delete[] WORKD;
//        delete[] WORKL;
//        delete[] RWORK;
        return 0;
    }
//
//    delete[] WORKD;
//    delete[] WORKL;
//    delete[] RWORK;

    return nconv;
}

int eigs_cn(const CMatType& A, CVecType& vals, CMatType& vecs,  int nev, std::string whch, double tol, int maxit, int ncv)
{
    size_t m=A.n_rows;
    assert(m==A.n_cols);

    auto MultAx=[&A,m](Complex in[], Complex out[]) -> void
    {
//        DOUT("before multiplication");
        CVecType invec(in,m,false), outvec(out,m,false);
        outvec = A*invec;
//        DOUT("after multiplication");
    };
    return eigs_cn(MultAx,m,vals,vecs,nev,whch,tol,maxit,ncv);
}

int eigs_rn(std::function<void (double*,double*)> MultOPx, int N, CVecType& vals, CMatType& vecs, int nev, std::string whch, double tol, int maxit, int ncv)
{
/// TODO (valentin#1#09/23/2013): replace all dynamically allocated arrays with memptr() from arma matrices/vectors

    vals.reset();
    vecs.reset();
    /// PARAMS FOR DNAUPD ---------------------------------------------------------------------------------------------------------//
    int mode=1; /// standard EV problem A*x = lam*x
    maxit=std::max(maxit,N);

    int IDO=0;
    char BMAT='I';
    char WHICH[3];
    strcpy(WHICH,whch.c_str());
    int NEV=nev;
    double TOL=tol;
    int NCV = (ncv==0) ? std::min(std::max(50,2*NEV+1),N) : ncv; /// from MATLAB:eigs:872
    int LDV = N;
    int LWORKL=3*NCV*(NCV + 2);
    int IPARAM[11] = {1,0,maxit,1,0,0,mode,0,0,0,0};
    int IPNTR[14];
    int INFONAUP=1;
//    int INFONAUP=0;

    RVecType RESID(N,fill::randn);
//    RVecType RESID(N);
    RMatType V(LDV,NCV);
    RVecType WORKD(3*N);
    RVecType WORKL(LWORKL);
//    double * WORKD = new double[3*N];
//    double * WORKL = new double[LWORKL];
//    uint ct=0;
    /// DNAUPD -------------------------------------------------------------------------------------------------------------------------//
    while (IDO!=99)
    {
        dnaupd_(&IDO,&BMAT,&N,WHICH,&NEV,&TOL,RESID.memptr(),&NCV,V.memptr(),&LDV,IPARAM,IPNTR,WORKD.memptr(),WORKL.memptr(),&LWORKL,&INFONAUP);

        switch (IDO)
        {
        case -1:
            cerr<<"-1 not implemented"<<endl;
            break;
        case 1: /// compute Z = B * X and Y = OP * Z
            /// For Simple EV Problem Y = Z
            MultOPx(&WORKD[IPNTR[0]-1],&WORKD[IPNTR[1]-1]);
            break;
        case 2:
            cerr<<"2 not implemented"<<endl;
            break;
        case 3:
            cerr<<"3 not implemented"<<endl;
            break;
        case 4:
            cerr<<"4 not implemented"<<endl;
            break;
        case 99:/// ARPACK HAS CONVERGED
            break;
        default:
            cerr<<"IDO has unknown value"<<endl;
        }
    }
    int nconv=0;

    if (IDO==99 && INFONAUP==0)
    {
        nconv=IPARAM[4];
    }
    else
    {
        cerr<<"no convergence in DNAUPD, on exit: "<<INFONAUP<<endl;

//        delete[] WORKD;
//        delete[] WORKL;
        return 0;
    }


    /// PARAMS FOR DNEUPD -------------------------------------------------------------------------------------------------------------------------//
    int INFONEUP=0;
    int RVEC=1;
    char HOWMNY='A';
    IVecType SELECT(NCV,fill::zeros);
//    SELECT.zeros();
    RVecType EVR(NEV+1);
    RVecType EVI(NEV+1);
    RMatType Zmat(N,NEV+1);
//    double * DR = EVR.memptr();
//    double * DI = EVI.memptr();
//    double * Z = Zmat.memptr();
    int LDZ=N;
    double SIGMAR=0., SIGMAI=0.;
    RVecType WORKEV(3*NCV);


    /// DNEUPD -------------------------------------------------------------------------------------------------------------------------------------//

//    dneupd_(&RVEC,&HOWMNY,SELECT.memptr(),DR,DI,Z,&LDZ,&SIGMAR,&SIGMAI,WORKEV.memptr(),&BMAT,&N,WHICH,&NEV,&TOL,RESID.memptr(),&NCV,V.memptr(),&LDV,IPARAM,IPNTR,WORKD,WORKL,&LWORKL,&INFONEUP);
    dneupd_(&RVEC,&HOWMNY,SELECT.memptr(),EVR.memptr(),EVI.memptr(),Zmat.memptr(),&LDZ,&SIGMAR,&SIGMAI,WORKEV.memptr(),
            &BMAT,&N,WHICH,&NEV,&TOL,RESID.memptr(),&NCV,V.memptr(),&LDV,IPARAM,IPNTR,WORKD.memptr(),WORKL.memptr(),&LWORKL,&INFONEUP);

    if(INFONEUP==0)
    {
        nconv=IPARAM[4];
        size_t newsize=nconv;

        EVR.resize(newsize);
        EVI.resize(newsize);

        vals=CVecType(EVR,EVI); /// fill in eigenvalues

        if (RVEC==1)
        {
            vecs.set_size(N,newsize);
            size_t i=0;
            while (i<newsize)
            {
                if (std::abs(EVI(i))<1e-14) /// real eigenvalue
                {
                    vecs.col(i)=CVecType(Zmat.col(i),zeros(N));
                    ++i;
                }
                else /// complex pair
                {
                    vecs.col(i)=CVecType(Zmat.col(i),Zmat.col(i+1));
                    vecs.col(i+1)=CVecType(Zmat.col(i),-Zmat.col(i+1));
                    i+=2;
                }
            }
        }
        else vals=sort(vals,1);
    }
    else
    {
        cerr<<"no convergence in DNEUP, on exit: "<<INFONEUP<<endl;

//        delete[] WORKD;
//        delete[] WORKL;
        return 0;
    }

    /// CLEAN UP
//    delete[] WORKD;
//    delete[] WORKL;

    return nconv;
}


//template<typename VT>
int eigs_rn(const RMatType& A, CVecType& vals, CMatType& vecs,  int nev, std::string whch, double tol, int maxit, int ncv)
{
    size_t m=A.n_rows;
    assert(m==A.n_cols);

    auto MultAx=[&A,m](double in[], double out[])
    {
        RVecType invec(in,m,false), outvec(out,m,false);
        outvec = A*invec;
    };

    return eigs_rn(MultAx,m,vals,vecs,nev,whch,tol,maxit,ncv);
}





