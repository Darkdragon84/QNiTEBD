#ifndef EIGS_H_
#define EIGS_H_

#include <functional>
#include <iostream>
#include <assert.h>

#include "arma_typedefs.h"

extern "C"
{
//    struct /// variables controlling debugging in ARPACK, see debuc.doc in DOCUMENTS for detail
//    {
//        int logfil, ndigit, mgetv0=0, msaupd=0, msaup2=0, msaitr=0, mseigt=0, msapps=0, msgets=0, mseupd=0, mnaupd=1, mnaup2=0, mnaitr=0, mneigh=0, mnapps=0, mngets=3, mneupd=4, mcaupd=0, mcaup2=0, mcaitr=0, mceigh=0, mcapps=0, mcgets=0, mceupd=0;
//    } DEBUG_;

    void dnaupd_(int* IDO, char* BMAT, int* N, char WHICH[], int* NEV, double* TOL, double RESID[], int* NCV, double V[], int* LDV, int IPARAM[],
                 int IPNTR[], double WORKD[], double WORKL[], int* LWORKL, int* INFO);
    void dneupd_(int* RVEC, char* HOWMNY, int SELECT[], double DR[], double DI[], double Z[], int* LDZ, double* SIGMAR, double* SIGMAI, double WORKEV[], char* BMAT, int* N,
                 char WHICH[], int* NEV, double* TOL, double RESID[], int* NCV, double V[], int* LDV, int IPARAM[], int IPNTR[], double WORKD[], double WORKL[], int* LWORKL, int* INFO);
    void znaupd_(int* IDO, char* BMAT, int* N, char WHICH[], int* NEV, double* TOL, Complex RESID[], int* NCV, Complex V[], int* LDV, int IPARAM[],
                 int IPNTR[], Complex WORKD[], Complex WORKL[], int* LWORKL, double RWORK[], int* INFO);
    void zneupd_(int* RVEC, char* HOWMNY, int SELECT[], Complex D[], Complex Z[], int* LDZ, Complex* SIGMA, Complex WORKEV[], char* BMAT, int* N, char WHICH[], int* NEV,
                 double* TOL, Complex RESID[], int* NCV, Complex V[], int* LDV, int IPARAM[], int IPNTR[], Complex WORKD[], Complex WORKL[], int* LWORKL, double RWORK[], int* INFO);
}

int eigs_rn(std::function<void (double*,double*)> MultOPx, int N, CVecType& vals, CMatType& vecs, int nev, std::string whch="LM", double tol=1e-14, int maxit=500, int ncv=0);
int eigs_rn(const RMatType& A, CVecType& vals, CMatType& vecs,  int nev, std::string whch="LM", double tol=1e-14, int maxit=500, int ncv=0);

int eigs_cn(std::function<void (Complex*,Complex*)> MultOPx, int n, CVecType& vals, CMatType& vecs, int nev, std::string whch="LM", double tol=1e-14, int maxit=500, int ncv=0);
int eigs_cn(const CMatType& A, CVecType& vals, CMatType& vecs,  int nev, std::string whch="LM", double tol=1e-14, int maxit=500, int ncv=0);

//template<typename VT>
//void ShowArray(VT arr[], size_t len){for (size_t i=0; i<len; ++i)cout<<arr[i]<<endl;}
#endif // EIGS_H_
