#ifndef __LONGRANGESPINLADDER_H
#define __LONGRANGESPINLADDER_H
#include "hams.h"
#include "fitting.h"
#include "hambuilder.h"

class LongRangeSpinLadder : public MPOBuilder
    {
    public:

    LongRangeSpinLadder(const Model& model_,
                        const ExpFit& fit1, const ExpFit& fitXY, const ExpFit& fitZ)
        : 
        MPOBuilder(model_),
        fit1_(fit1),
        fitXY_(fitXY),
        fitZ_(fitZ),
        f_(0)
        { }

    LongRangeSpinLadder(const Model& model_, const Callable& f,
                        const ExpFit& fit1, const ExpFit& fitXY, const ExpFit& fitZ)
        : 
        MPOBuilder(model_),
        fit1_(fit1),
        fitXY_(fitXY),
        fitZ_(fitZ),
        f_(&f)
        { }

    operator const IQMPO&() { init(); return QH; }

    operator const MPO&() 
        { 
        init(); 
        if(H.isNull())
            H = QH.toMPO();
        return H; 
        }

    private:

    /////////////
    //
    // Data Members

    const ExpFit &fit1_,
                 &fitXY_,
                 &fitZ_;
    IQMPO QH;
    MPO H;

    const Callable* f_;

    //
    /////////////

    void
    init()
        {
        QH = IQMPO(model);
        const int N = model.NN();

        //Determine number of complex chi's
        int ncomplex1 = 0; 
        std::vector<bool> iscomplex1(fit1_.nexp()+1,false);
        for(int k = 1; k <= fit1_.nexp(); ++k) 
        if(fabs(fit1_.ImChi()(k)) > 1E-15)
            {
            iscomplex1.at(k) = true;
            ncomplex1 += 1;
            }
        std::cout << "ncomplex1 = " << ncomplex1 << std::endl;

        int ncomplexXY = 0; 
        std::vector<bool> iscomplexXY(fitXY_.nexp()+1,false);
        for(int k = 1; k <= fitXY_.nexp(); ++k) 
        if(fabs(fitXY_.ImChi()(k)) > 1E-15)
            {
            iscomplexXY.at(k) = true;
            ncomplexXY += 1;
            }
        std::cout << "ncomplexXY = " << ncomplexXY << std::endl;

        int ncomplexZ = 0; 
        std::vector<bool> iscomplexZ(fitZ_.nexp()+1,false);
        for(int k = 1; k <= fitZ_.nexp(); ++k) 
        if(fabs(fitZ_.ImChi()(k)) > 1E-15)
            {
            iscomplexZ.at(k) = true;
            ncomplexZ += 1;
            }
        std::cout << "ncomplexZ = " << ncomplexZ << std::endl;

        const int p1 = fit1_.nexp()+ncomplex1,
                  pXY = fitXY_.nexp()+ncomplexXY,
                  pZ = fitZ_.nexp()+ncomplexZ,
                  ko1 = 2*p1,
                  koXY = 2*pXY,
                  koZ = 2*pZ,
                  ds = 2*ko1+2*koXY+1,    //start of diagonal part
                  kd = ko1+koZ+2,    //bond dimension of diagonal part
                  k  = (ds-1)+kd;    //total bond dimension

        //std::cout << "p_ = " << p_ << std::endl;
        //std::cout << "ncomplex = " << ncomplex << std::endl;
        //std::cout << "ko = " << ko << std::endl;
        //std::cout << "ds = " << ds << std::endl;
        //std::cout << "kd = " << kd << std::endl;
        //std::cout << "k  = " << k << std::endl;

        std::vector<IQIndex> iqlinks(N+1);

        //The names of these indices refer to their Nf quantum numbers (plus or minus), 
        //but they can have various sz quantum numbers depending on the type of site they follow
        std::vector<Index> q0(N+1),
                           qP(N+1), 
                           qM(N+1);

        for(int i = 0; i <= N; ++i)
            {
            qP.at(i) = Index(nameint("qP_",i),ko1+koXY);
            qM.at(i) = Index(nameint("qM_",i),ko1+koXY);
            q0.at(i) = Index(nameint("q0_",i),kd);

            iqlinks.at(i) = IQIndex(nameint("hl",i),
                                    qP[i],QN(-2),
                                    qM[i],QN(+2),
                                    q0[i],QN( 0));
            }

        std::vector<Index> start_inds(1); 
        start_inds[0] = q0.at(0);

        std::vector<Index> end_inds(1); 
        end_inds[0] = q0.at(N);

        for(int j = 1; j <= N; ++j)
            {
            //Create j^th A (an IQTensor)
            IQTensor &W = QH.AAnc(j);
            IQIndex row = conj(iqlinks.at(j-1)),
                    col = iqlinks.at(j);

            int leg = (j%2==1 ? 1 : 2);

            const Real xj = (j/2)+1;

            W = IQTensor(conj(model.si(j)),model.siP(j),row,col);

            //Identity string operators
            W += model.id(j) * row(ds) * col(ds);
            W += model.id(j) * row(k) * col(k);

            //Long range Heisenberg interactions along legs
            for(int type = 1; type <= 3; ++type)
                {
                int r = 0;
                IQTensor start_op, end_op;
                if(type == 1)
                    {
                    //S+ S- interactions
                    r = 1;
                    //std::cout << "Starting r = 2 = " << r << std::endl;
                    start_op = model.sp(j)*0.5;
                    end_op = model.sm(j);
                    }
                else if(type == 2)
                    {
                    r = ko1+1;
                    //std::cout << "Starting r = ko1+2 = " << r << std::endl;
                    start_op = model.sm(j)*0.5;
                    end_op = model.sp(j);
                    }
                else if(type == 3)
                    {
                    r = ds+1;
                    //std::cout << "Starting r = ds+2 = " << r << std::endl;
                    start_op = model.sz(j);
                    end_op = model.sz(j);
                    }

                if(f_ != 0)
                    {
                    start_op *= (*f_)(xj);
                    end_op *= (*f_)(xj);
                    }

                //Iterate over legs: a = 1,2
                for(int a = 1; a <= 2; ++a)
                    {
                    bool this_leg = (a == leg);
                    //Iterate over exponentials in the fit l = 1,2,...,p
                    for(int l = 1; l <= fit1_.nexp(); ++l) 
                        {
                        if(this_leg)
                            {
                            W += model.id(j)*row(r)*col(r);
                            W += end_op*row(r)*col(ds);
                            W += fit1_.ReChi()(l)*start_op*row(k)*col(r);
                            }
                        else
                            {
                            W += fit1_.ReLambda()(l)*model.id(j)*row(r)*col(r);
                            }
                        ++r;
                        if(iscomplex1.at(l))
                            {
                            if(this_leg)
                                {
                                W += model.id(j)*row(r)*col(r);
                                W += end_op*row(r)*col(ds);
                                W += -fit1_.ImChi()(l)*start_op*row(k)*col(r);
                                }
                            else
                                {
                                W += fit1_.ReLambda()(l)*model.id(j)*row(r)*col(r);
                                W += -fit1_.ImLambda()(l)*model.id(j)*row(r-1)*col(r);
                                W += fit1_.ImLambda()(l)*model.id(j)*row(r)*col(r-1);
                                }
                            ++r;
                            }
                        } //for l
                    } //for a

                } //for type

            //Interactions between legs
            for(int type = 1; type <= 3; ++type)
                {
                int r = 0;
                IQTensor start_op, end_op;
                const ExpFit* pfit = 0;
                const std::vector<bool>* piscomplex = 0;

                if(type == 1)
                    {
                    //S+ S- interactions
                    r = 2*ko1+1;
                    //std::cout << "Starting r = 2 = " << r << std::endl;
                    start_op = model.sp(j)*0.5;
                    end_op = model.sm(j);
                    pfit = &fitXY_;
                    piscomplex = &iscomplexXY;
                    }
                else if(type == 2)
                    {
                    r = 2*ko1+koXY+1;
                    //std::cout << "Starting r = ko1+2 = " << r << std::endl;
                    start_op = model.sm(j)*0.5;
                    end_op = model.sp(j);
                    pfit = &fitXY_;
                    piscomplex = &iscomplexXY;
                    }
                else if(type == 3)
                    {
                    r = ds+ko1+1;
                    //std::cout << "Starting r = ds+2 = " << r << std::endl;
                    start_op = model.sz(j);
                    end_op = model.sz(j);
                    pfit = &fitZ_;
                    piscomplex = &iscomplexZ;
                    }

                if(f_ != 0)
                    {
                    start_op *= (*f_)(xj);
                    end_op *= (*f_)(xj);
                    }

                //Iterate over legs: a = 1,2
                for(int a = 1; a <= 2; ++a)
                    {
                    bool this_leg = (a == leg);
                    //Iterate over exponentials in the fit l = 1,2,...,p
                    for(int l = 1; l <= pfit->nexp(); ++l) 
                        {
                        if(this_leg)
                            {
                            W += pfit->ReLambda()(l)*model.id(j)*row(r)*col(r);
                            if(leg == 1)
                                W += pfit->ReChi()(l)*start_op*row(k)*col(r);
                            else
                                W += pfit->ReLambda()(l)*pfit->ReChi()(l)*start_op*row(k)*col(r);
                            }
                        else
                            {
                            W += model.id(j)*row(r)*col(r);
                            W += pfit->ReLambda()(l)*end_op*row(r)*col(ds);
                            }
                        ++r;
                        if(piscomplex->at(l))
                            {
                            if(this_leg)
                                {
                                W += pfit->ReLambda()(l)*model.id(j)*row(r)*col(r);
                                W += -pfit->ImLambda()(l)*model.id(j)*row(r-1)*col(r);
                                W += pfit->ImLambda()(l)*model.id(j)*row(r)*col(r-1);

                                if(leg == 1)
                                    {
                                    W += -pfit->ImChi()(l)*start_op*row(k)*col(r);
                                    }
                                else
                                    {
                                    W += -pfit->ImLambda()(l)*pfit->ImChi()(l)*start_op*row(k)*col(r-1);
                                    W += -pfit->ImLambda()(l)*pfit->ReChi()(l)*start_op*row(k)*col(r);
                                    W += -pfit->ReLambda()(l)*pfit->ImChi()(l)*start_op*row(k)*col(r);
                                    }
                                }
                            else
                                {
                                W += model.id(j)*row(r)*col(r);
                                W += -pfit->ImLambda()(l)*end_op*row(r-1)*col(ds);
                                W += (pfit->ImLambda()(l)+pfit->ReLambda()(l))*end_op*row(r)*col(ds);
                                }
                            ++r;
                            }
                        } //for l
                    } //for a

                } //for type

            }

        QH.AAnc(1) = makeLedge(iqlinks.at(0),start_inds) * QH.AA(1);
        QH.AAnc(N) = QH.AA(N) * makeRedge(conj(iqlinks.at(N)),end_inds); 

        }

    };

void inline
HTerms(const SpinHalf& model,
       Real LambdaXY, Real LambdaZ, std::vector<IQMPO>& terms,
       const Option& opt1 = Option(), const Option& opt2 = Option())
    {
    const int N = model.NN();
    const int nx = N/2;

    OptionSet oset(opt1,opt2);

    int cutoff = oset.intOrDefault("Cutoff",-1);
    if(oset.defined("Cutoff"))
        std::cout << boost::format("Using cutoff %d for Hamiltonian\n") % cutoff << std::endl;
    else
        std::cout << "Not using a cutoff for the Hamiltonian" << std::endl;

    HamBuilder hb(model);

    IQMPO term; //to be used as a temporary throughout

    terms.clear();

    std::cout << boost::format("terms.size() = %d") % terms.size() << std::endl;

    if(oset.defined("Offset"))
        {
        terms.push_back(term);
        std::cout << "Pushed back null term" << std::endl;
        }

    int count = 0;

    //Power-law interaction on chains
    std::cout << "Creating intra-chain interactions" << std::endl;
    for(int c = 1; c <= 2; ++c)
    for(int r1 = 1; r1 <= nx; ++r1)
        {
        //std::cout << boost::format("r1 = %d/%d\n") % r1 % nx << std::endl;
        for(int r2 = r1+1; r2 <= nx; ++r2)
            {
            std::vector<IQMPO> subterms;
            //std::cout << boost::format("PL int %d of %d\n") % r1 % nx << std::endl;
            
            if(cutoff >= 0 && (r2-r1) > (cutoff+1)) continue;

            int i = 0, j = 0;
            if(c == 1)
                {
                i = 2*r1-1; 
                j = 2*r2-1;
                }
            else
                {
                i = 2*r1; 
                j = 2*r2;
                }


            std::cout << boost::format("i = %d, j = %d, term # %d") % i % j % (++count) << std::endl;

            Real fac = 1./pow(r2-r1,3);

            /*
            if(r1 == 1)
                {
                Real dist = fabs(1.*(r2-r1));
                std::cout << boost::format("Intra-chain: r1 = %d, r2 = %d, dist = %d, fac = %.3f\n") % r1 % r2 % dist % fac << std::endl;
                }
                */

            hb.getMPO(i,model.sz(i),j,model.sz(j),
                      term,fac);
            subterms.push_back(term);
            hb.getMPO(i,model.sp(i),j,model.sm(j),
                      term,0.5*fac);
            subterms.push_back(term);
            hb.getMPO(i,model.sm(i),j,model.sp(j),
                      term,0.5*fac);
            subterms.push_back(term);

            sum(subterms,term);
            terms.push_back(term);
            //std::cout << boost::format("Now terms.size() = %d") % terms.size() << std::endl;
            }
        }

    //Inter chain interactions
    std::cout << "Creating inter-chain interactions" << std::endl;
    for(int c = 1; c <= 2; ++c)
    for(int r1 = 1; r1 <= nx; ++r1)
        {
        //std::cout << boost::format("r1 = %d/%d\n") % r1 % nx << std::endl;
        for(int r2 = r1; r2 <= nx; ++r2)
            {
            std::vector<IQMPO> subterms;

            if(c == 2 && r2 == r1) continue;

            if(cutoff >= 0 && (r2-r1) > cutoff) continue;

            int i = 0, j = 0;
            if(c == 1)
                {
                i = 2*r1-1; 
                j = 2*r2;
                }
            else
                {
                i = 2*r1; 
                j = 2*r2-1;
                }
            std::cout << boost::format("i = %d, j = %d, term # %d") % i % j % (++count) << std::endl;

            Real d = fabs(1.*(r2-r1));

            Real lxy = 1./pow(d*d+1,1.5)*(1-(1.-LambdaXY)/(d*d+1));
            Real lz = 1./pow(d*d+1,1.5)*(1.-(1.-LambdaZ)/(d*d+1));

            /*
            if(r1 == 1)
                {
                Real dist = fabs(1.*(r2-r1));
                std::cout << boost::format("Inter-chain: r1 = %d, r2 = %d, dist = %d, lxy = %.3f, lz = %.3f\n") 
                             % r1 % r2 % dist % lxy % lz << std::endl;
                }
                */


            if(lz != 0)
                {
                hb.getMPO(i,model.sz(i),j,model.sz(j),
                          term,lz);
                subterms.push_back(term);
                }

            if(lxy != 0)
                {
                hb.getMPO(i,model.sp(i),j,model.sm(j),
                          term,0.5*lxy);
                subterms.push_back(term);
                hb.getMPO(i,model.sm(i),j,model.sp(j),
                          term,0.5*lxy);
                subterms.push_back(term);
                }

            if(!subterms.empty())
                {
                sum(subterms,term);
                terms.push_back(term);
                //std::cout << boost::format("Now terms.size() = %d") % terms.size() << std::endl;
                }

            }
        }
    }


#endif
