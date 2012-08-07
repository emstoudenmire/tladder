#ifndef __NNSPINLADDER_H
#define __NNSPINLADDER_H
#include "hams.h"
#include "fitting.h"
#include "hambuilder.h"

class NNSpinLadder : public MPOBuilder
    {
    public:

    NNSpinLadder(const Model& model_,
                 Real LambdaXY, Real LambdaZ)
        : 
        MPOBuilder(model_),
        LambdaZ_(LambdaZ),
        LambdaXY_(LambdaXY)
        { }

    operator const IQMPO&() { init(); return H; }

    private:

    /////////////
    //
    // Data Members

    Real LambdaZ_,
         LambdaXY_;

    IQMPO H;

    //
    /////////////

    void
    init()
        {
        H = IQMPO(model);
        const int N = model.NN();


        const int kk = 2,
                  ds = 5,    //start of diagonal part
                  kd = kk+2, //bond dimension of diagonal part
                  k  = (ds-1)+kd;    //total bond dimension

        std::vector<IQIndex> iqlinks(N+1);

        //The names of these indices refer to their Nf quantum numbers (plus or minus), 
        //but they can have various sz quantum numbers depending on the type of site they follow
        std::vector<Index> q0(N+1),
                           qP(N+1), 
                           qM(N+1);

        for(int i = 0; i <= N; ++i)
            {
            qP.at(i) = Index(nameint("qP_",i),kk);
            qM.at(i) = Index(nameint("qM_",i),kk);
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
            IQTensor &W = H.AAnc(j);
            IQIndex row = conj(iqlinks.at(j-1)),
                    col = iqlinks.at(j);

            int leg = (j%2==1 ? 1 : 2);

            W = IQTensor(conj(model.si(j)),model.siP(j),row,col);

            //Identity string operators
            W += model.id(j) * row(ds) * col(ds);
            W += model.id(j) * row(k) * col(k);

            //S+ S- terms
            if(leg == 1)
                {
                W += model.sp(j) * row(k) * col(1) * (LambdaXY_/2.);
                }
            W += model.sp(j) * row(k) * col(2) * 0.5;
            W += model.id(j) * row(2) * col(1);
            W += model.sm(j) * row(1) * col(ds);

            //S- S+ terms
            if(leg == 1)
                {
                W += model.sm(j) * row(k) * col(3) * (LambdaXY_/2.);
                }
            W += model.sm(j) * row(k) * col(4) * 0.5;
            W += model.id(j) * row(4) * col(3);
            W += model.sp(j) * row(3) * col(ds);

            //Sz Sz terms
            if(leg == 1)
                {
                W += model.sz(j) * row(k) * col(6) * LambdaZ_;
                }
            W += model.sz(j) * row(k) * col(7);
            W += model.id(j) * row(7) * col(6);
            W += model.sz(j) * row(6) * col(ds);
            }

        H.AAnc(1) = makeLedge(iqlinks.at(0),start_inds) * H.AA(1);
        H.AAnc(N) = H.AA(N) * makeRedge(conj(iqlinks.at(N)),end_inds); 
        }

    };


#endif
