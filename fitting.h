#ifndef __FITTING_H
#define __FITTING_H

#include "itensor.h"

Vector inline
fabs(const Vector& v) 
    { 
    Vector res(v); 
    for(int i = 1; i <= v.Length(); ++i) 
        res(i) = fabs(res(i)); 
    return res; 
    }

Real inline
max(const Vector& v) 
    { 
    Real res = -1E10; 
    for(int i = 1; i <= v.Length(); ++i) 
        res = max(res,v(i)); 
    return res; 
    }


void inline 
PseudoInverse(const Matrix& M, Matrix& Minv, Real cutoff = 0)
    {
    Vector D; Matrix U,V;
    SVD(M,U,D,V,1E-2);

    //{
    //Matrix DD(D.Length(),D.Length());
    //DD.Diagonal() = D;
    //Matrix err = M - (U*DD*V);
    //Real sumerrsq = Trace(err * err.t());
    //std::cout << boost::format("SVD err is %.2E\n") % sqrt(sumerrsq/(M.Nrows()*M.Ncols()));
    //}

    int m = D.Length();
    Matrix Lambda(m,m); Lambda = 0.0;
    for(int i = 1; i <= m; ++i)
        {
        if(fabs(D(i)/D(1)) > cutoff) Lambda(i,i) = 1./D(i);
        }

    Minv = V.t()*Lambda*U.t();
    }

class Cplx
    {
public:
    Real a,b;
    Cplx(Real a_, Real b_) : a(a_), b(b_) { }

    void 
    operator*=(const Cplx& rhs)
        {
        Real na = a*rhs.a - b*rhs.b;
        b = b*rhs.a + a*rhs.b;
        a = na;
        }

    Cplx operator*(const Cplx& rhs)
        { 
        Cplx res(*this); 
        res *= rhs; 
        return res; 
        }

    void operator+=(const Cplx& rhs) { a += rhs.a; b += rhs.b; }

    Cplx pow(int d)
        {
        if(d < 0) Error("Cplx::pow only accepts positive powers.");
        Cplx res(*this);
        for(int n = 1; n <= d-1; ++n) res *= (*this);
        return res;
        }
    };


class Callable
    {
    public:

    Callable() { }

    Real 
    operator()(Real x) const { return call(x); }

    private:

    Real virtual
    call(Real x) const = 0;

    };



class ExpFit
    {
    public:

    ExpFit();

    ExpFit(const Callable& f, int N, int nexp,
           const Option& opt1 = Option(), const Option& opt2 = Option());

    Real
    operator()(int d) const;

    void
    stats(bool& lambda_is_real, Real& maxdiff, 
          int& maxdiff_pos, Real& avgdiff) const;

    const Vector&
    ReLambda() const { return ReLambda_; }
    const Vector&
    ImLambda() const { return ImLambda_; }

    const Vector&
    ReChi() const { return ReChi_; }
    const Vector&
    ImChi() const { return ImChi_; }

    int
    nexp() const { return nexp_; }

    private:

    ///////////////////
    //
    // Data Members
    
    const Callable* f_;
    
    int N_,
        Nb_,
        nexp_;

    Vector ReLambda_,
           ImLambda_,
           ReChi_,
           ImChi_;

    //
    //////////////////

    void
    init();

    };

inline ExpFit::
ExpFit()
    :
    f_(0),
    N_(0),
    Nb_(0),
    nexp_(0)
    { }

inline ExpFit::
ExpFit(const Callable& f,int N, int nexp,
       const Option& opt1, const Option& opt2)
    :
    f_(&f),
    N_(N),
    Nb_(N-1),
    nexp_(nexp)
    { 
    OptionSet oset(opt1,opt2);

    //
    // Do automatic fit if requested
    //
    if(oset.defined("Auto"))
        {
        bool lambda_is_real = true;
        Real maxdiff = 100000.0;
        int maxdiff_pos = -1;
        Real avgdiff = 0;

        Real min_maxdiff = 10000.0;
        int best_nexp = 2;
        for(nexp_ = 2; nexp_ <= min(Nb_-2,nexp); ++nexp_)
            {
            init();
            stats(lambda_is_real,maxdiff,maxdiff_pos,avgdiff);
            if(!(oset.boolOrDefault("Quiet",false)))
                std::cout << "Trying nexp = " << nexp_ << ", maxdiff = " << maxdiff << std::endl;
            if(fabs(maxdiff) < min_maxdiff)
                {
                min_maxdiff = maxdiff;
                best_nexp = nexp_;
                }
            }

        nexp_ = best_nexp;
        }

    init();
    }


void inline ExpFit::
init()
    {
    if(nexp_ <= 0)
        {
        nexp_ = 0;
        ReChi_.ReDimension(1);
        ImChi_.ReDimension(1);
        ReChi_(1) = 0.0;
        ImChi_(1) = 0.0;
        ReLambda_.ReDimension(1);
        ImLambda_.ReDimension(1);
        ReLambda_(1) = 1.0;
        ImLambda_(1) = 0.0;
        return;
        }
    if(nexp_ > Nb_) error("Number of fitting parameters can't exceed number of sites.");

    Matrix M(Nb_-nexp_+1,nexp_);
    M = 0;

    for(int i = 0; i < nexp_; ++i)
    for(int j = 1; j <= Nb_-nexp_+1; ++j)
        {
        M(j,i+1) = (*f_)(i+j);
        }

    int m = M.Nrows();

    Matrix Q,R;

    int method = 3;
    if(method == 1)
        {
        Matrix Q1,R1;
        QRDecomp(M,Q1,R1);
        Orthog(Q1,nexp_,2);
        Matrix M1 = Q1.t()*M;

        Matrix Q2;
        QRDecomp(M1,Q2,R);

        Q = Q1*Q2;
        }
    else if(method == 2)
        {
        QRDecomp(M,Q,R);
        }
    else if(method == 3)
        {
        Vector D;
        Matrix V;
        SVD(M,Q,D,V);
        }

    //{
    //Matrix err = M-Q*R;
    //Real sumerrsq = Trace(err * err.t());
    //std::cout << boost::format("QR err is %.2E\n") % sqrt(sumerrsq/(p*m));
    //}

    //
    // Q1 is all but the last row of Q
    // Q2 all but the first row of Q
    //
    Matrix Q1 = Q.SubMatrix(1,m-1,1,nexp_),
           Q2 = Q.SubMatrix(2,m,1,nexp_);

    Matrix Q1i;
    PseudoInverse(Q1,Q1i);

    Matrix Meff = Q1i*Q2;

    GenEigenValues(Meff,ReLambda_,ImLambda_);

    //
    // Do least squares fit to
    // compute weights Chi
    //

    //Set up least squares problem:
    // L * chi = fv
    Matrix L(Nb_,2*nexp_); L = 0;
    Vector fv(Nb_); fv = 0;
    for(int i = 1; i <= Nb_; ++i)
        {
        fv(i) = (*f_)(i);
        for(int j = 1; j <= nexp_; ++j)
            {
            Cplx z = (Cplx(ReLambda_(j),ImLambda_(j))).pow(i);
            L(i,2*j-1) = z.a;
            L(i,2*j) = -z.b;
            }
        }

    Matrix Li;
    PseudoInverse(L,Li,1E-12);

    Vector ChiCombined = Li*fv;

    //Unpack Chi vector into real
    //and imaginary parts
    ReChi_.ReDimension(nexp_);
    ImChi_.ReDimension(nexp_);
    for(int j = 1; j <= nexp_; ++j)
        {
        ReChi_(j) = ChiCombined(2*j-1);
        ImChi_(j) = ChiCombined(2*j);
        }

    }

Real inline ExpFit::
operator()(int d) const
    {
    Cplx res(0,0);
    for(int k = 1; k <= ReLambda_.Length(); ++k) 
        {
        Cplx chi(ReChi_(k),ImChi_(k));
        Cplx lambda(ReLambda_(k),ImLambda_(k));
        res += chi*lambda.pow(d);
        }
    if(fabs(res.b) > 1E-3) 
        {
        static int count = 0;
        if(++count < 10)
            {
            std::cout << "WARNING: non-zero imaginary part of fit" << std::endl;
            std::cout << boost::format("    res.b = %.3E\n") % res.b;
            }
        }
    return res.a;
    }

void inline ExpFit::
stats(bool& lambda_is_real, Real& maxdiff, int& maxdiff_pos, Real& avgdiff) const
    {
    if(Norm(ImLambda_) > 1E-10) lambda_is_real = false;
    else lambda_is_real = true;

    avgdiff = 0.0;
    maxdiff = -1.0;
    maxdiff_pos = -1;
    Real mindiff = 10000.0;
    int mindiff_pos = -1;
    for(int d = 1; d <= Nb_; ++d)
        {
        Real diff = fabs((*f_)(d) - operator()(d));
        if(diff > maxdiff)
            {
            maxdiff = diff;
            maxdiff_pos = d;
            }
        if(diff < mindiff) 
            {
            mindiff = diff;
            mindiff_pos = d;
            }
        avgdiff += diff;
        }
    avgdiff /= Nb_;
    }


#endif
