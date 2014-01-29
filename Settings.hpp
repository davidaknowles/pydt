#define _USE_MATH_DEFINES

#ifndef SETTINGS_H
#define SETTINGS_H

#include <iostream>
#include <list>
#include <cstdlib>
#include <math.h>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/poisson_distribution.hpp>
#include <boost/random/discrete_distribution.hpp>
#include <string>
//#include <limits>

using namespace std;

class BaseSettings {
public:
  bool debug; 
  
  double theta, alpha;
    int D;

    int PoissonEventsOnFinalEdge;

    bool ResampleAttachmentLocation;
    bool ResampleDetachmentLocation; // TODO implement this
    bool sampleTheta, sampleAlpha, sampleC;
    int NumSamplesForImputation;

    double subtreeInflationFactor; // increase this to increase the average size of detached subtrees

    double c;
    double starFactor;
    boost::random::mt19937 gen;

    double nodeWidth;

    BaseSettings() {
      debug=false; 
        theta = 1.0;
        alpha = 0.0;
        sampleTheta = true;
        sampleAlpha = true;
        sampleC = true;
        PoissonEventsOnFinalEdge = 3;
        NumSamplesForImputation=3;
        subtreeInflationFactor = 2.0;
        ResampleAttachmentLocation = true;
        ResampleDetachmentLocation = false;
        //gen.seed(static_cast<unsigned int>(time(0)));
        gen.seed(0);
        c = 1.0;
        starFactor = 1.1;
        nodeWidth = 0.1;
    }

    // TODO: handle simpler cases efficiently
    double logDivergenceRateFactor(int m, int numInSubtree=1)
    {
        double res= 0;
        for (int i=0; i<numInSubtree; i++)
            res += -lgamma((double)m+(double)i-alpha)+lgamma((double)m+(double)i+1.0+theta);
        return res;
    }

    bool multifurcating()
    {
        return (theta+2.0*alpha) != 0.0;
    }

    // Integral of a
    virtual double A(double t) = 0;

    // Integral of a*
    virtual double Astar(double t) = 0;

    // Inverse of A
    virtual double invA(double x) = 0;

    // Inverse of A*
    virtual double invAstar(double x) = 0;

    // Inverse of A(t)+A*(t)
    virtual double invAplusAstar(double x) = 0;

    // Rate function
    virtual double a(double t) = 0;

    // Upper bound on rate function
    virtual double astar(double t) = 0;

    // sample uniform between fMin and fMax
    double fRand(double fMin = 0.0, double fMax = 1.0)
    {
        boost::uniform_real<> dist(fMin, fMax);
        return dist(gen);
    }

    double singleTerm(int i)
    {
        return lgamma((double)i-alpha)-lgamma((double)i+1.0+theta);
    }

  std::map<int,double> memo; 

    // Generalised Harmonic number
    double H_internal(int n) // memoize?
    {
        double result = 0;
        for (int i=1; i<=n; i++)
            result += exp(singleTerm(i));
        return result;
    }

  double H(int n){
    std::map<int,double>::iterator f = memo.find(n); 
    if (f == memo.end()){
      double newH=H_internal(n); 
      memo[n]=newH; 
      return newH; 
    } else {
      return *f; 
    }
  }

  void HypersChanged(){
    memo.clear(); 
  }

    // Sample an integer between min (inclusive) and max (inclusive)
    int iRand(int min, int max)
    {
        boost::random::uniform_int_distribution<> dist(min, max);
        return dist(gen);
    }

    // Sample a Poisson rv with mean
    int rpoiss(double mean)
    {
        if (mean==0.0) return 0;
        boost::random::poisson_distribution<> dist(mean);
        return dist(gen);
    }

    // Sample from a discrete distribution. Input must be normalised
    int rdiscrete(boost::numeric::ublas::vector<double> a)
    {
        boost::uniform_real<> dist(0.0,1.0);
        double u= dist(gen);
        for (int i=0; i<a.size(); i++)
        {
            u-=a[i];
            if (u<0.0)
                return i;
        }
        throw 1;
    }

    // Sample from hist
    int sampleMap(map<int,int> &hist, int total)
    {
        int u = iRand(0, total);
        for (map<int,int>::iterator i = hist.begin(); i != hist.end(); i++)
        {
            u -= i->second;
            if (u < 0)
                return i->first;
        }
        throw 1;
    }

};

class Settings : public BaseSettings
{
public:

    Settings() : BaseSettings()
    {

    }

    virtual double A(double t)
    {
        return -c*log(1.0-t);
    };

    // must have Astar(t) > A(t)_
    virtual double Astar(double t)
    {
        return A(t)*starFactor;
    }

    virtual double invA(double x)
    {
        return 1.0-exp(-x/c);
    };

    virtual double invAstar(double x)
    {
        return invA(x/starFactor);
    }

    virtual double invAplusAstar(double x)
    {
        return invA(x/(starFactor-1.0));
    }

    virtual double a(double t)
    {
        return c / (1.0 - t);
    }

    virtual double astar(double t)
    {
        return a(t)*starFactor;
    }
};

class BoundedSettings : public BaseSettings
{
public:

    BoundedSettings() : BaseSettings()
    {

    }

    virtual double A(double t)
    {
        return c*t;
    };

    // must have Astar(t) > A(t)_
    virtual double Astar(double t)
    {
        return A(t)*starFactor;
    }

    virtual double invA(double x)
    {
        return x/c;
    };

    virtual double invAstar(double x)
    {
        return invA(x/starFactor);
    }

    virtual double invAplusAstar(double x)
    {
        return invA(x/(starFactor-1.0));
    }

    virtual double a(double t)
    {
        return c;
    }

    virtual double astar(double t)
    {
        return a(t)*starFactor;
    }
};
#endif
