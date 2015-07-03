#ifndef TYPEDEFS_ARM_H_
#define TYPEDEFS_ARM_H_

#include <assert.h>
#include <complex>
#include <map>
#include <memory>
#include <armadillo>
#include <time.h>

#ifdef SYMS
#define _USE_SYMMETRIES_
#endif // SYMS

#ifdef VERBOSE
#define DOUT(str) do { std::cout << str << std::endl; } while( false )
#else
#define DOUT(str) do { } while ( false )
#endif

using namespace arma;

/// typedefs
typedef double Real;
typedef unsigned int uint;
typedef std::complex<double> Complex;

typedef Mat<Real> RMatType;
typedef Mat<Complex> CMatType;
typedef Col<Real> RVecType;
typedef Col<Complex> CVecType;
typedef Col<int> IVecType;
typedef Col<uint> UIVecType;

typedef std::shared_ptr<RMatType> pRMatType;
typedef std::shared_ptr<CMatType> pCMatType;
typedef std::shared_ptr<RVecType> pRVecType;
typedef std::shared_ptr<CVecType> pCVecType;
typedef std::shared_ptr<IVecType> pIVecType;

#ifdef COMPLEX
using Scalar = Complex;
using VecType = CVecType;
using MatType = CMatType;
#else
using Scalar = Real;
using VecType = RVecType;
using MatType = RMatType;
#endif

/// additional typedefs
enum dirtype {r,l,s,c};
typedef std::pair<uint,uint> uiipair;
typedef std::map<uint,uint> uiimap;

#endif // TYPEDEFS_ARM_H_
