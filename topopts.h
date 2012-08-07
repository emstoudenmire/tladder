#ifndef __TOPOPTS_H
#define __TOPOPTS_H
#include "DMRGObserver.h"

#define Format boost::format
#define Cout std::cout
#define Endl std::endl

class TopOpts : public DMRGObserver
    {
    public:

    typedef DMRGObserver 
    Parent;

    TopOpts(const IQMPS& psi, const SpinHalf& model, const std::string& pfix = "");

    const std::string& 
    prefix() const { return prefix_; }
    void 
    prefix(std::string val) { prefix_ = val; }

    int 
    notifyTimes() const { return notify_times_; }
    void 
    notifyTimes(int val) { notify_times_ = val; }

    Real 
    esAccuracy() const { return es_accuracy_; }
    void 
    esAccuracy(Real val) { es_accuracy_ = val; }

    virtual void 
    measure(int sw, int ha, int b, const SVDWorker& svd, Real energy,
            const Option& opt1 = Option(), const Option& opt2 = Option(),
            const Option& opt3 = Option(), const Option& opt4 = Option());

    bool virtual
    checkDone(int sw, const SVDWorker& svd, Real energy,
              const Option& opt1 = Option(), const Option& opt2 = Option());

    private:

    /////////////
    //
    // Data Members
    //

    const IQMPS& psi_;
    const SpinHalf& model_;

    bool do_plot_self;
    int notify_times_;
       
    std::string prefix_;

    Real last_es_,
         curr_es_,
         es_accuracy_;

    //
    /////////////

    };

inline TopOpts::
TopOpts(const IQMPS& psi, const SpinHalf& model, const std::string& pfix)
    : 
    psi_(psi), 
    model_(model), 
    do_plot_self(true),
    notify_times_(-1), 
    prefix_(pfix),
    last_es_(1000),
    curr_es_(-1000),
    es_accuracy_(-1)
    { }

void inline TopOpts::
measure(int sw, int ha, int b, const SVDWorker& svd, Real energy,
        const Option& opt1, const Option& opt2,
        const Option& opt3, const Option& opt4)
    {
    Parent::measure(sw,ha,b,svd,energy);

    if(ha == 2)
    if(abs(b - model_.NN()/2) <= 1)
        {
        const std::string pstring = (prefix_ == "" ? "" : prefix_ + "_");

        Vector eigs = svd.eigsKept(model_.NN()/2);

        Real splitting = 0;
        for(int j = 1; j <= eigs.Length()/2; ++j)
            splitting += sqrt(fabs(eigs(2*j-1)))-sqrt(fabs(eigs(2*j)));

        std::cout << boost::format("\nEntanglement splitting at bond %d = %.10f") % b % splitting << std::endl;

        if(b == model_.NN()/2)
            {
            last_es_ = curr_es_;
            curr_es_ = splitting;
            }
        }


    //Notify where the DMRG worker is every 'notify_times' sites
    if(notify_times_ > 0 && (b+1)%(psi_.NN()/notify_times_) == 0)
        std::cout << "In sweep " << sw << ", reached bond " << b << std::endl;


    /*
    if(b == 3 && ha == 2)
        {
        std::string wfname_ = "gs_psi";
        std::cout << "Writing wavefunction to disk (as " << wfname_ << ")" << std::endl;
        writeToFile(wfname_,psi_);
        }
    */

    }

bool TopOpts::
checkDone(int sw, const SVDWorker& svd, Real energy,
          const Option& opt1, const Option& opt2)
    {
    if(Parent::checkDone(sw,svd,energy,opt1,opt2)) return true;

    if(es_accuracy_ > 0 && fabs(last_es_-curr_es_) < es_accuracy_)
        {
        Cout << Format("\nEntanglement splitting accuracy goal met, %.2E < %.2E.\n") 
                % (fabs(last_es_-curr_es_)) % es_accuracy_ << Endl;
        return true;
        }

    return false;
    }

#undef Format
#undef Cout
#undef Endl

#endif
