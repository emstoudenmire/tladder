#ifndef __PARAMS_H
#define __PARAMS_H
#include "Sweeps.h"
#include "DMRGObserver.h"
#include "input.h"

class Params
    {
    public:

    Real
    cutoff,
    esaccuracy,
    J,
    K,
    LambdaXY,
    LambdaZ,
    orth_weight,
    param_end,
    param_start,
    param_step,
    pinning,
    xi;

    int
    do_param_sweep,
    do_plot_self,
    do_timing,
    interaction_cutoff,
    max_p,
    max_p_leg,
    max_p_rung,
    maxm,
    minm,
    min_sweeps,
    nn,
    nstates,
    nsweeps,
    nwarm,
    nx,
    p,
    printH,
    quiet,
    quiet_dmrg,
    smooth,
    stagger_pinning,
    triplet_sector,
    use_tmpdir,
    write_m;

    std::string
    nthreads,
    runmode,
    sweep_param,
    sweep_scheme,
    wfname,
    write_dir;

    Sweeps::Scheme scheme;

    Params() {}

    Params(InputGroup& basic) 
        {
        //Get mandatory params
        basic.GetIntM("nx",nx);

        //Defaults for optional params

        //Real
        cutoff = 1E-8;
        esaccuracy = -1;
        J = 1;
        LambdaXY = 1;
        LambdaZ = 1;
        orth_weight = 1;
        param_end = -37;
        param_start = -37;
        param_step = -1;
        pinning = 0;
        xi = 1;

        //int
        do_param_sweep = 0;
        do_plot_self = 0;
        do_timing = 0;
        interaction_cutoff = -1;
        max_p = 25;
        max_p_leg = -1;
        max_p_rung = -1;
        maxm = 1000;
        minm = 1;
        nn = 0;
        nstates = 1;
        nsweeps = 5;
        nwarm = nsweeps-1;
        min_sweeps = 2;
        p = -1;
        printH = 0;
        quiet = 1;
        quiet_dmrg = 1;
        smooth = 0;
        stagger_pinning = 0;
        use_tmpdir = 0;
        write_m = -1;

        //string
        nthreads = "";
        runmode = "solve";
        sweep_param = "lambdaxy";
        sweep_scheme = "ramp_m";
        wfname = "";
        write_dir = "";

        //Get optional params
        basic.GetYesNo("do_plot_self",do_plot_self);
        basic.GetYesNo("do_timing",do_timing);
        basic.GetReal("esaccuracy",esaccuracy);
        basic.GetReal("J",J);
        basic.GetReal("K",K);
        basic.GetReal("LambdaXY",LambdaXY);
        basic.GetReal("LambdaZ",LambdaZ);
        basic.GetInt("interaction_cutoff",interaction_cutoff);
        basic.GetInt("max_p",max_p);
        basic.GetInt("max_p_leg",max_p_leg);
        basic.GetInt("max_p_rung",max_p_rung);
        basic.GetInt("min_sweeps",min_sweeps);
        basic.GetYesNo("nn",nn);
        basic.GetInt("nstates",nstates);
        basic.GetInt("nsweeps",nsweeps);
        basic.GetString("nthreads",nthreads);
        basic.GetInt("nwarm",nwarm);
        basic.GetReal("orth_weight",orth_weight);
        basic.GetInt("p",p);
        basic.GetReal("param_end",param_end);
        basic.GetReal("param_start",param_start);
        basic.GetReal("param_step",param_step);
        basic.GetReal("pinning",pinning);
        basic.GetYesNo("printH",printH);
        basic.GetYesNo("quiet",quiet);
        basic.GetYesNo("quiet_dmrg",quiet_dmrg);
        basic.GetString("runmode",runmode);
        basic.GetYesNo("smooth",smooth);
        basic.GetYesNo("stagger_pinning",stagger_pinning);
        basic.GetString("sweep_scheme",sweep_scheme);
        basic.GetYesNo("triplet_sector",triplet_sector);
        basic.GetYesNo("use_tmpdir",use_tmpdir);
        basic.GetString("wfname",wfname);
        basic.GetString("write_dir",write_dir);
        basic.GetInt("write_m",write_m);
        basic.GetReal("xi",xi);

        //Derived parameters
        scheme = Sweeps::ramp_m;
        if(sweep_scheme == "ramp_m") scheme = Sweeps::ramp_m;
        else if(sweep_scheme == "fixed_cutoff") scheme = Sweeps::fixed_cutoff;
        else if(sweep_scheme == "fixed_m") scheme = Sweeps::fixed_m;
        else if(sweep_scheme == "exp_m") scheme = Sweeps::exp_m;
        else if(sweep_scheme == "table") scheme = Sweeps::table;
        else Error("Sweep scheme unrecognized.");

        if(param_start != param_end || param_step != -1)
            do_param_sweep = 1;


        if(write_m > 1)
            Global::options().add(WriteM(write_m));

        if(use_tmpdir)
            {
            const char* tempdir = getenv("TMPDIR");
            std::cout << "Setting write directory to " << tempdir << std::endl;
            Global::options().add(WriteDir(tempdir));
            }
        else
        if(write_dir != "")
            {
            std::cout << "Setting write directory to " << write_dir << std::endl;
            Global::options().add(WriteDir(write_dir));
            }

        if(nthreads != "")
            {
            std::cout << "Setting thread env variables to " + nthreads << std::endl;
            const char* ntcs = nthreads.c_str();

            setenv("VECLIB_MAXIMUM_THREADS",ntcs,1);
            char* vnt = getenv("VECLIB_MAXIMUM_THREADS");
            if(vnt != NULL)
                std::cout << "VECLIB_MAXIMUM_THREADS = " << vnt << std::endl;
            else
                std::cout << "VECLIB_MAXIMUM_THREADS not defined" << std::endl;

            setenv("OMP_NUM_THREADS",ntcs,1);
            const char* ont = getenv("OMP_NUM_THREADS");
            if(ont != NULL)
                std::cout << "OMP_NUM_THREADS = " << ont << std::endl;
            else
                std::cout << "OMP_NUM_THREADS not defined" << std::endl;
            }
        }

    ~Params() {}
    };

#endif
