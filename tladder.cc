#include "params.h"
#include "core.h"
#include "model/spinhalf.h"
#include "model/spinone.h"
#include "hams/heisenberg.h"
#include "writedata.h"
#include "fitting.h"
#include "LongRangeSpinLadder.h"
#include "NNSpinLadder.h"
#include "topopts.h"
using boost::format;
using namespace std;

Params params;

bool 
fexist(const string& filename) 
    { 
	ifstream file(filename.c_str());
	return file.good();
    }
bool 
fexist(const format& ffname) 
    { 
    return fexist(ffname.str());
    }

class Dipole : public Callable
    {
    public:

    Dipole() { }

    private: 

    Real virtual
    call(Real d) const
        {
        return 1./pow(d,3);
        }

    };

class InterLeg : public Callable
    {
    public:

    InterLeg(Real Lambda) 
        :
        Lambda_(Lambda)
        { }

    private: 

    Real Lambda_;

    Real virtual
    call(Real d) const
        {
        Real x = (d-1);
        return 1./pow(x*x+1,1.5)*(1-(1-Lambda_)/(x*x+1));
        }

    };

class TanhSmoothing : public Callable
    {
    public:

    TanhSmoothing(Real nx, Real xi) 
        :
        nx_(nx),
        xi_(fabs(xi))
        { 
        if(xi_ == 0)
            Error("Can't set xi to zero");
        }

    private:

    Real nx_,
         xi_;

    Real
    call(Real j) const
        {
        Real dj = 2+(j < nx_/2 ? (j-1) : nx_-j);
        return y(1-dj/xi_);
        }

    Real
    y(Real x) const
        {
        if(x <= 0)
            return 1;
        else if(x >= 1)
            return 0;
        else
            return 0.5*(1-tanh((x-0.5)/(x*(1-x))));
        }

    };

void
makeLongRangeH(const Model& model, IQMPO& H)
    {
    if(params.nn)
        Error("NN Hamiltonian requested");

    const int nx = params.nx;
    const int N = 2*nx; //2 leg ladder
    const Real LambdaXY = params.LambdaXY;
    const Real LambdaZ = params.LambdaZ;

    TanhSmoothing smoothing(nx,params.xi);

    //
    // Compute long range fits
    //
    ExpFit fit,fitXY,fitZ;

    int p = params.p;
    int max_p_leg = (params.max_p_leg == -1 ? params.max_p : params.max_p_leg);
    int max_p_rung = (params.max_p_rung == -1 ? params.max_p : params.max_p_rung);

    Dipole f;
    fit = ExpFit(f,nx,(p < 0 ? max_p_leg : p),Auto(p < 0),Quiet());

    InterLeg lxy(LambdaXY);
    fitXY = ExpFit(lxy,nx,(p < 0 ? max_p_rung : p),Auto(p < 0),Quiet());

    InterLeg lz(LambdaZ);
    fitZ = ExpFit(lz,nx,(p < 0 ? max_p_rung : p),Auto(p < 0),Quiet());
    Real totZ1 = 0, totZ2 = 0;
    for(int n = 1; n <= fitZ.ReChi().Length(); ++n)
        {
        totZ1 += fitZ.ReChi()(n);
        totZ2 += fitZ.ReChi()(n)*fitZ.ReLambda()(n);
        }
    //cout << format("LambdaZ = %.10f, totZ1 = %.10f, totZ2 = %.10f\n") % LambdaZ % totZ1 % totZ2 << endl;

    Real avgdiff = 0.0;
    Real maxdiff = -1.0;
    int maxdiff_pos = -1;
    bool lambda_is_real = true;

    fit.stats(lambda_is_real,maxdiff,maxdiff_pos,avgdiff);
    cout << "Dipole fit:" << endl;
    cout << "    Number of exponentials = " << fit.nexp() << endl;
    cout << "    maxdiff = " << format("%.3E") % maxdiff << ", at " << maxdiff_pos << endl;
    cout << "    avgdiff = " << format("%.3E") % avgdiff << endl;
    cout << endl;

    fitXY.stats(lambda_is_real,maxdiff,maxdiff_pos,avgdiff);
    cout << "XY fit:" << endl;
    cout << "    Number of exponentials = " << fitXY.nexp() << endl;
    cout << "    maxdiff = " << format("%.3E") % maxdiff << ", at " << maxdiff_pos << endl;
    cout << "    avgdiff = " << format("%.3E") % avgdiff << endl;
    cout << endl;

    fitZ.stats(lambda_is_real,maxdiff,maxdiff_pos,avgdiff);
    cout << "Z fit:" << endl;
    cout << "    Number of exponentials = " << fitZ.nexp() << endl;
    cout << "    maxdiff = " << format("%.3E") % maxdiff << ", at " << maxdiff_pos << endl;
    cout << "    avgdiff = " << format("%.3E") % avgdiff << endl;
    cout << endl;

    Vector fitted_V2(nx), exact_V2(nx);
        {
        for(int j = 1; j <= nx; ++j) 
            {
            fitted_V2(j) = fit(j);
            exact_V2(j) = f(j);
            }
        writedata("fitted_V2",fitted_V2,1,params.do_plot_self);
        writedata("exact_V2",exact_V2,1,params.do_plot_self);
        Vector V2_diff = exact_V2-fitted_V2;
        writedata("V2_diff",V2_diff,1,params.do_plot_self);
        }

        {
        for(int j = 1; j <= nx; ++j) 
            {
            fitted_V2(j) = fitXY(j);
            exact_V2(j) = lxy(j);
            }
        writedata("XY_fitted_V2",fitted_V2,1,params.do_plot_self);
        writedata("XY_exact_V2",exact_V2,1,params.do_plot_self);
        Vector V2_diff = exact_V2-fitted_V2;
        writedata("XY_V2_diff",V2_diff,1,params.do_plot_self);
        }

        {
        for(int j = 1; j <= nx; ++j) 
            {
            fitted_V2(j) = fitZ(j);
            exact_V2(j) = lz(j);
            }
        writedata("Z_fitted_V2",fitted_V2,1,params.do_plot_self);
        writedata("Z_exact_V2",exact_V2,1,params.do_plot_self);
        Vector V2_diff = exact_V2-fitted_V2;
        writedata("Z_V2_diff",V2_diff,1,params.do_plot_self);
        }

    if(params.smooth)
        {
        cout << "\nUsing smooth long range model.\n" << endl;
        H = LongRangeSpinLadder(model,smoothing,fit,fitXY,fitZ);
        }
    else
        {
        cout << "\nUsing long range model.\n" << endl;
        H = LongRangeSpinLadder(model,fit,fitXY,fitZ);
        }

    }

int main(int argc, char* argv[])
    {
    //Get parameter file
    if(argc != 2)
        {
        cout << "Usage: " << argv[0] << " inputfile." << endl;
        return 0;
        }
    string infilename(argv[1]);
    InputFile infile(infilename);
    InputGroup basic(infile,"basic");
    basic.quiet = true;

    params = Params(basic);

    const int nx = params.nx;
    const int N = 2*nx; //2 leg ladder
    const Real LambdaXY = params.LambdaXY;
    const Real LambdaZ = params.LambdaZ;


    string model_name = (format("model_%d") % nx).str();

    SpinHalf model;
    if(fexist(model_name))
        {
        cout << "Reading model " << model_name << " from disk." << endl;
        readFromFile(model_name,model);
        }
    else
        {
        model = SpinHalf(N);
        }

    TanhSmoothing smoothing(nx,params.xi);
    if(params.smooth)
        {
        Vector F(nx);
        for(int j = 1; j <= nx; ++j)
            F(j) = smoothing(j);
        writedata("F",F,1,params.do_plot_self);
        }



    InputGroup table(basic,"sweeps");
    Sweeps sweeps(params.nsweeps,table);
    cout << sweeps;

    if(params.runmode == "solve")
    {

    IQMPS psi(model);
    if(params.wfname != "" && fexist(params.wfname))
        {
        cout << "Reading wavefunction " << params.wfname << " from file." << endl;
        readFromFile(params.wfname,psi);
        }
    else
        {
        InitState initState(N);
        for(int i = 1; i <= N; ++i) 
            initState(i) = (i%2==1 ? model.Dn(i) : model.Up(i));
        cout << "Creating new initial state wavefunction." << endl;
        psi = IQMPS(model,initState);
        }

    Real En = 0;
    /*
    if(params.interaction_cutoff > -1)
        {
        cout << "params.interaction_cutoff = " << params.interaction_cutoff << endl;
        vector<IQMPO> terms;
        if(1)
            {
            HTerms(model,LambdaXY,LambdaZ,terms,Cutoff(params.interaction_cutoff));
            cout << format("Summing %d terms to form H") % terms.size() << endl;
            IQMPO Hs;
            sum(terms,Hs);
            cout << "Starting DMRG" << endl;
            En = dmrg(psi,Hs,sweeps,opts,Quiet(params.quiet_dmrg));
            cout << format("GS Energy = %.10f\n") % En;

            Real checkEn = 0;
            Foreach(const IQMPO& term, terms)
                {
                checkEn += psiHphi(psi,term,psi);
                }
            cout << format("Check Energy = %.10f\n") % checkEn;

            cout << format("Energy using full H = %.10f\n") % psiHphi(psi,H,psi);
            }
        else
            {
            HTerms(model,LambdaXY,LambdaZ,terms,Cutoff(params.interaction_cutoff),Offset());
            cout << format("Starting DMRG using %d H terms") % (terms.size()-1) << endl;
            En = dmrg(psi,terms,sweeps,opts,Quiet(params.quiet_dmrg));
            }
        }
        */



    IQMPO H;

    Real *swp_paramP = 0;
    if(params.sweep_param == "lambdaxy")
        swp_paramP = &(params.LambdaXY);

    Real &swp_param = *swp_paramP;

    if(!params.do_param_sweep)
        {
        params.param_start = swp_param;
        params.param_end = swp_param;
        params.param_step = 1000;
        }

    for(swp_param = params.param_start; (swp_param-params.param_end) < 1E-8; swp_param += params.param_step)
        {
        cout << format("\n\nMaking Hamiltonian with sweep param %s = %.10f\n") % params.sweep_param % swp_param << endl;
        if(params.nn)
            {
            cout << "\nUsing nearest-neighbor model.\n" << endl;
            H = NNSpinLadder(model,params.LambdaXY,params.LambdaZ);
            }
        else
            {
            makeLongRangeH(model,H);
            }

        TopOpts opts(psi,model);
        if(params.esaccuracy > 0)
            opts.esAccuracy(params.esaccuracy);
        //opts.notifyTimes(4);

        En = dmrg(psi,H,sweeps,opts,Quiet(params.quiet_dmrg));
        cout << format("GS Energy = %.10f\n") % En;

        writeToFile(format("gs_psi_%s_%.4f")%params.sweep_param%swp_param,psi);
        writeToFile(model_name,model);
        }



    /*
    vector<IQMPO> terms;
    HTerms(model,LambdaXY,LambdaZ,terms,Cutoff(params.interaction_cutoff));

    Real checkEn = 0;
    //int tot = terms.size();
    int count = 1;
    Foreach(const IQMPO& term, terms)
        {
        //cout << format("Term %d of %d") % count % tot << endl;
        checkEn += psiHphi(psi,term,psi);
        ++count;
        }
    cout << format("Check Energy = %.10f\n") % checkEn;
    */

    } //end runmode solve
    else
    if(params.runmode == "gap")
    {
    const int nstates = params.nstates;

    IQMPO H;
    if(params.nn)
        {
        cout << "\nUsing nearest-neighbor model.\n" << endl;
        H = NNSpinLadder(model,params.LambdaXY,params.LambdaZ);
        }
    else
        {
        makeLongRangeH(model,H);
        }

    vector<IQMPS> psi;
    psi.reserve(nstates);
    vector<Real> energy(nstates+1,1000);

    //Make initial state (Neel state)
    InitState initState(N);
    for(int i = 1; i <= N; ++i) 
        initState(i) = (i%2==1 ? model.Dn(i) : model.Up(i));

    Matrix olap(nstates,nstates);
    olap = 1;
    Matrix Heff(nstates,nstates);
    Heff = 0;

    Vector eigs;
    Matrix U;

    for(int state = 0; state < nstates; ++state)
        {
        IQMPS newpsi(model,initState);

        TopOpts opts(newpsi,model);

        cout << format("\n\nBeginning DMRG calculation for state %d\n") % state << endl;

        if(state == 0)
            {
            energy.at(state) = dmrg(newpsi,H,sweeps,opts,Quiet(params.quiet_dmrg));
            }
        else
            {
            energy.at(state) = dmrg(newpsi,H,psi,sweeps,opts,
                                    Weight(params.orth_weight),Quiet(params.quiet_dmrg));
            }

        cout << "Energies:" << endl;
        for(int s = 0; s <= state; ++s)
            {
            cout << format("   %d  %.10f\n") % s % energy.at(s);
            }

        cout << "Updating overlap matrix..." << endl;
        for(int s = 0; s < state; ++s)
            {
            olap.el(s,state) = psiphi(psi.at(s),newpsi);
            olap.el(state,s) = olap.el(s,state);
            }

        cout << "Updating Heff matrix..." << endl;
        for(int s = 0; s < state; ++s)
            {
            Heff.el(s,state) = psiHphi(psi.at(s),H,newpsi);
            Heff.el(state,s) = Heff.el(s,state);
            }
        Heff.el(state,state) = energy.at(state);

        cout << "Current overlap matrix:" << endl;
        for(int r = 1; r <= state+1; ++r)
            {
            for(int c = 1; c <= state+1; ++c)
                {
                if(r == c)
                    cout << "1.000000 ";
                else
                    cout << format("%.2E ") % fabs(olap(r,c));
                }
            cout << endl;
            }

        cout << endl;

        cout << "Current Heff matrix:" << endl;
        for(int r = 1; r <= state+1; ++r)
            {
            for(int c = 1; c <= state+1; ++c)
                {
                cout << format("%+.2E ") % Heff(r,c);
                }
            cout << endl;
            }

        cout << "Eigenvalues of Heff matrix:" << endl;
        EigenValues(Heff,eigs,U);
        for(int s = 0; s < nstates; ++s)
            {
            cout << format("   %d  %.10f\n") % s % eigs.el(s);
            }

        Matrix HT = U.t()*Heff*U;
        cout << "Transformed Heff matrix:" << endl;
        for(int r = 1; r <= state+1; ++r)
            {
            for(int c = 1; c <= state+1; ++c)
                {
                cout << format("%+.2E ") % HT(r,c);
                }
            cout << endl;
            }

        Matrix olapT = U.t()*olap*U;
        cout << "Transformed overlap matrix:" << endl;
        for(int r = 1; r <= state+1; ++r)
            {
            for(int c = 1; c <= state+1; ++c)
                {
                cout << format("%+.2E ") % olapT(r,c);
                }
            cout << endl;
            }

        psi.push_back(newpsi);
        }

    } //end runmode gap
    else
    {
    Error("Runmode not recognized");
    }

    return 0;
    }
