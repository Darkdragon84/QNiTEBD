#ifndef ITERATIVE_H
#define ITERATIVE_H

#include <iostream>
#include <functional>

#include "arma_typedefs.h"


template<typename T>
uint gmres(const std::function<Col<T> (const Col<T>&)>& Afun, const Col<T>& b, Col<T>& x, const Col<T>& x0, Real tol, uint maxit, uint krylovdim=0, bool verbose=false)
{
    uint flag = 0;
    uint dim = b.n_elem;
    assert(x0.n_elem == dim);
    if (maxit > dim)
    {
        if (verbose) cout<<"maxit > dim="<<dim<<", setting back to maxit="<<dim<<endl;
        maxit = dim;
    }
    if (krylovdim==0) krylovdim = min(dim,maxit);
    if (krylovdim > dim)
    {
        if (verbose) cout<<"krylovdim > dim="<<dim<<", setting back to krylovdim="<<dim<<endl;
        krylovdim = dim;
    }

    Col<T> r0 = b - Afun(x0);
    x = x0;
    Col<T> r(r0), w;
    Real bnorm = norm(b), rnorm = norm(r), thresh=1e-15;
    T beta,htmp1,htmp2,phase;
    Mat<T> V(dim,krylovdim+1);
    Mat<T> H(krylovdim+1,krylovdim+1);
    Col<T> c(krylovdim+1),s(krylovdim+1), gamma(krylovdim+1);
    gamma(0) = rnorm;
    V.col(0) = r/rnorm;

//    ct = 0;
    uint kdim=0;
    for (uint j = 0;j < krylovdim; ++j)
    {
//        cout<<j+1<<":"<<endl;
        /// modified Gram-Schmidt (use qr?)
        w = Afun(V.col(j));
        for (uint i=0;i<j+1;++i)
        {
            H(i,j) = cdot(V.col(i),w);
            w -= H(i,j)*V.col(i);
        }
        H(j+1,j)=norm(w);
        /// apply (all j previous) Givens rotations in current column j to update H such that H(j+1,j) can be made zero below
        for (uint i=0;i<j;++i)
        {
            htmp1 = H(i,j);
            htmp2 = H(i+1,j);
            H(i,j)   =        c(i)*htmp1 + s(i)*htmp2;
            H(i+1,j) = -conj(s(i))*htmp1 + c(i)*htmp2;
        }
        beta = norm(H(span(j,j+1),j));
//        beta = sqrt(H(j,j)*H(j,j) + H(j+1,j)*H(j+1,j)); /// that basically means that w didn't change after subtracting the last v, i.e. w would be linearly dependent on last v.
        if (abs(beta) < thresh) break;
        else
        {
            V.col(j+1) = w/H(j+1,j);
            /** calculate new Givens parameters to make H(j+1,j)=0
                we want to find parameters c and s, with c real and s complex, such that
                |  c   s | | h(j  ,j) |  = | r |
                | -s*  c | | h(j+1,j) |  = | 0 |
                and c^2 + |s|^2 = 0, i.e. the above is a unitary with the special choice c real.
                The below choice does it with r = |h(j,j)|/h(j,j) * beta
            **/
            phase = H(j,j)/abs(H(j,j));
            s(j) = phase*conj(H(j+1,j))/beta;
            c(j) = abs(H(j,j))/beta;
//            cout<<H(j+1,j)*c(j) - H(j,j)*conj(s(j))<<" should be zero"<<endl;
            /// we are overwriting H with R column after column. In the last iteration R is an upper triangular matrix with all zeros on the bottom
            H(j,j) = phase*beta;
            H(j+1,j) = 0;
            gamma(j+1) = -conj(s(j))*gamma(j);
            gamma(j) *= c(j); /// important to do that after setting gamma(j+1)!!
            ++kdim;
        }
        rnorm = abs(gamma(j+1));
//        if(verbose) cout<<j+1<<":"<<rnorm<<endl;
        if (rnorm < tol*bnorm)
        {
            flag = kdim;
            break;
        }
    }
    /// after completion, reconstruct x from alpha and v
    /// determine alpha from solving R*alpha=gamma. Since R is upper triangular this can be done efficiently from the bottom up

    if (H.n_cols > kdim) H.resize(kdim,kdim);
    if (gamma.n_elem > kdim) gamma.resize(kdim);
    if (V.n_cols > kdim) V.resize(dim,kdim);

//    cout<<"V is ONS:"<<norm(V.t()*V - eye<RMatType>(V.n_cols,V.n_cols),"fro")<<endl;

    Col<T> alpha = solve(trimatu(H),gamma);
    x += V*alpha;

    if (verbose)
    {
        r = b - Afun(x);
        Real relnorm = norm(r)/bnorm;
        if (flag > 0) cout<<"converged after "<<flag<<" steps to |r|/|b|="<<relnorm<<endl;
        else cerr<<"no convergence after "<<kdim<<" steps, |r|/|b|="<<relnorm<<endl;
    }
    return flag;
}

#endif
