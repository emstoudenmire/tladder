#ifndef __WRITEDATA_H
#define __WRITEDATA_H
#include <string>
#include <iostream>
#include <fstream>
#include "matrix.h"
#include "boost/format.hpp"

inline void System(const char* cstr) { int res = system(cstr); res++; }
inline void System(const std::string str) { System(str.c_str()); }
inline void System(const std::stringstream& s) { System(s.str()); }
inline void System(const boost::format fmt) { System(fmt.str()); }

//Vector data types; x, y and error data
void inline
writedata(const char* cstr, const Vector& xdat,const Vector& ydat, const Vector& edat, bool do_plot_self = false)
    {
    std::ofstream f(cstr);
    if(xdat.Length() != ydat.Length()) Error("xdat and ydat Lengths don't match.");
    for(int j = 1; j <= xdat.Length(); ++j) f << boost::format("%.20f %.20f %.20f\n") % xdat(j) % ydat(j) % edat(j);
    f.close();
    //if(do_plot_self) System(boost::format("$HOME/tools/plot_self %s") % cstr);
    if(do_plot_self) System(boost::format("plot_self %s") % cstr);
    }
inline void writedata(const std::string str, const Vector& xdat, const Vector& ydat, const Vector& edat, bool do_plot_self = false) { writedata(str.c_str(),xdat,ydat,edat,do_plot_self); }
inline void writedata(const std::stringstream& s,const Vector& xdat, const Vector& ydat, const Vector& edat, bool do_plot_self=false) { writedata(s.str(),xdat,ydat,edat,do_plot_self); }
inline void writedata(const boost::format fmt, const Vector& xdat, const Vector& ydat, const Vector& edat, bool do_plot_self=false) { writedata(fmt.str(),xdat,ydat,edat,do_plot_self); }

//Vector data types; x, y data
inline void writedata(const char* cstr, const Vector& xdat,const Vector& ydat, bool do_plot_self = false)
{
    Vector edat = ydat; edat = 0;
    writedata(cstr,xdat,ydat,edat,do_plot_self);
}
inline void writedata(const std::string str, const Vector& xdat, const Vector& ydat, bool do_plot_self = false) { writedata(str.c_str(),xdat,ydat,do_plot_self); }
inline void writedata(const std::stringstream& s,const Vector& xdat, const Vector& ydat, bool do_plot_self=false) { writedata(s.str(),xdat,ydat,do_plot_self); }
inline void writedata(const boost::format fmt, const Vector& xdat, const Vector& ydat, bool do_plot_self=false) { writedata(fmt.str(),xdat,ydat,do_plot_self); }

//std std::vector data types; x, y and error data
template<typename T>
void writedata(const char* cstr, const std::vector<T>& xdat,const std::vector<T>& ydat, const std::vector<T>& edat, bool do_plot_self = false)
{
    Vector X((int) xdat.size()); for(int n = 1; n <= X.Length(); ++n) X(n) = xdat.at(n-1);
    Vector Y((int) ydat.size()); for(int n = 1; n <= Y.Length(); ++n) Y(n) = ydat.at(n-1);
    Vector E((int) edat.size()); for(int n = 1; n <= E.Length(); ++n) E(n) = edat.at(n-1);
    writedata(cstr,X,Y,E,do_plot_self);
}
template<typename T>
inline void writedata(const std::string str, const std::vector<T>& xdat, const std::vector<T>& ydat, const std::vector<T>& edat, bool do_plot_self = false) { writedata(str.c_str(),xdat,ydat,edat,do_plot_self); }
template<typename T>
inline void writedata(const std::stringstream& s,const std::vector<T>& xdat, const std::vector<T>& ydat, const std::vector<T>& edat, bool do_plot_self=false) { writedata(s.str(),xdat,ydat,edat,do_plot_self); }
template<typename T>
inline void writedata(const boost::format fmt, const std::vector<T>& xdat, const std::vector<T>& ydat, const std::vector<T>& edat, bool do_plot_self=false) { writedata(fmt.str(),xdat,ydat,edat,do_plot_self); }

//std std::vector data types; x, y and error data
template<typename T>
void writedata(const char* cstr, const std::vector<T>& xdat,const std::vector<T>& ydat, bool do_plot_self = false)
{
    Vector X((int) xdat.size()); for(int n = 1; n <= X.Length(); ++n) X(n) = xdat.at(n-1);
    Vector Y((int) ydat.size()); for(int n = 1; n <= Y.Length(); ++n) Y(n) = ydat.at(n-1);
    writedata(cstr,X,Y,do_plot_self);
}
template<typename T>
inline void writedata(const std::string str, const std::vector<T>& xdat, const std::vector<T>& ydat, bool do_plot_self = false) { writedata(str.c_str(),xdat,ydat,do_plot_self); }
template<typename T>
inline void writedata(const std::stringstream& s,const std::vector<T>& xdat, const std::vector<T>& ydat, bool do_plot_self=false) { writedata(s.str(),xdat,ydat,do_plot_self); }
template<typename T>
inline void writedata(const boost::format fmt, const std::vector<T>& xdat, const std::vector<T>& ydat, bool do_plot_self=false) { writedata(fmt.str(),xdat,ydat,do_plot_self); }

template<typename T>
void writedata(const char* cstr,const std::vector<T>& ydat, bool do_plot_self = false)
{
    Vector X((int) ydat.size()); for(int n = 1; n <= X.Length(); ++n) X(n) = n-1;
    Vector Y((int) ydat.size()); for(int n = 1; n <= Y.Length(); ++n) Y(n) = ydat.at(n-1);
    writedata(cstr,X,Y,do_plot_self);
}
template<typename T>
inline void writedata(const std::string str, const std::vector<T>& ydat, bool do_plot_self = false) { writedata(str.c_str(),ydat,do_plot_self); }
template<typename T>
inline void writedata(const std::stringstream& s, const std::vector<T>& ydat, bool do_plot_self=false) { writedata(s.str(),ydat,do_plot_self); }
template<typename T>
inline void writedata(const boost::format fmt, const std::vector<T>& ydat, bool do_plot_self=false) { writedata(fmt.str(),ydat,do_plot_self); }

inline void writedata(const char* cstr, const Vector& dat, Real Delta = 1, bool do_plot_self = false)
{
    std::ofstream f(cstr);
    for(int j = 1; j <= dat.Length(); ++j) f << boost::format("%.20f %.20f\n") % (Delta*j) % dat(j);
    f.close();
    //if(do_plot_self) System(boost::format("$HOME/tools/plot_self %s") % cstr);
    if(do_plot_self) System(boost::format("plot_self %s") % cstr);
}
inline void writedata(const std::string str, const Vector& dat, Real Delta = 1, bool do_plot_self = false) { writedata(str.c_str(),dat,Delta,do_plot_self); }
inline void writedata(const std::stringstream& s,const Vector& dat, Real Delta = 1, bool do_plot_self=false) { writedata(s.str(),dat,Delta,do_plot_self); }
inline void writedata(const boost::format fmt, const Vector& dat, Real Delta = 1, bool do_plot_self=false) { writedata(fmt.str(),dat,Delta,do_plot_self); }

inline void writedata(const char* cstr, const Matrix& dat, bool do_plot_self = false)
{           
    std::ofstream f(cstr);
    for(int r = 1; r <= dat.Nrows(); ++r)
    for(int c = 1; c <= dat.Ncols(); ++c) 
        f << r SP c SP dat(r,c) << "\n";
    f.close();
    //if(do_plot_self) System(boost::format("$HOME/tools/plot_self %s") % cstr);
    if(do_plot_self) System(boost::format("plot_self %s") % cstr);
}           
inline void writedata(const std::string str, const Matrix& dat, bool do_plot_self=false) { writedata(str.c_str(),dat,do_plot_self); }
inline void writedata(const std::stringstream& s,const Matrix& dat, bool do_plot_self=false) { writedata(s.str(),dat,do_plot_self); }
inline void writedata(const boost::format fmt, const Matrix& dat, bool do_plot_self=false) { writedata(fmt.str(),dat,do_plot_self); }

inline void GetXYDataFromFile(Vector& xvalues, Vector& yvalues,  const char* filename   )
{
    std::vector<Real> xs,ys;
    std::ifstream infile(filename);
    std::cerr << "Reading file " << filename << " for data...\n\n" ;
    if(infile.is_open()) 
    {
        float x,y,lastx;
        lastx = 0;
        while(infile.good())
        {
            std::string line;
            getline(infile,line);
            sscanf(line.c_str(),"%f %f",&x,&y);
            std::cerr << boost::format("Got x=%f,y=%E\n") % x % y;
            
            //if (x != lastx) //for some reason, getline prints the last line twice
            xs.push_back((Real) x); ys.push_back((Real) y);
            //xvalues(i) = (Real) x;
            //yvalues(i) = (Real) y;
            
            lastx = x;
        }
        infile.close();
    }
    int Ns = ys.size() - 1;
	if (Ns != yvalues.Length()) Error("Data and system have different numbers of sites");
    std::cerr << "total # of lattice points read in = " <<  Ns  << "\n\n";
    std::cerr << "first x-value is " << xs.at(0) << " and last x value is " << xs.at(Ns-1) << "\n\n";
    //std::cerr << "first y-value is " << ys.at(0) << " and last y value is " << ys.at(Ns-1) << "\n\n";
    //std::cerr << "Ns value of x and y = " << xs.at(Ns) << " and " << ys.at(Ns) << "\n";
    //std::cerr << "Ns-2 value of x and y = "<< xs.at(Ns-2) << " and " << ys.at(Ns-2) << "\n";
    
    for(int j = 1; j <= Ns; ++j) {
        xvalues(j) = xs.at(j-1);
        yvalues(j) = ys.at(j-1);
    }
    std::cerr << "all finished" << std::endl;
    
}

inline void GetXYDataFromFile(Vector& xvalues, Vector& yvalues, const std::string filename)  {
    GetXYDataFromFile(xvalues, yvalues, filename.c_str());
}
inline void GetXYDataFromFile(Vector& xvalues, Vector& yvalues, const std::stringstream& filename)  {
    GetXYDataFromFile(xvalues, yvalues, filename.str());
}
inline void GetXYDataFromFile(Vector& xvalues, Vector& yvalues, const boost::format filename)  {
    GetXYDataFromFile(xvalues, yvalues, filename.str());
}

#endif
