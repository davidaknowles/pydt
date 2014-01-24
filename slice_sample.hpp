#define _USE_MATH_DEFINES

#ifndef SLICE_SAMPLE_H
#define SLICE_SAMPLE_H

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <math.h>

template<typename T>
double slice_sample(double xx, T& logdist, boost::random::mt19937 gen, bool step_out = true, double width=1.0)
{

    double log_Px = logdist(xx);

    boost::uniform_real<> dist(0.0,1.0);
    double log_uprime = log(dist(gen)) + log_Px;

    double xprime = xx;
    double rr = dist(gen);
    double x_l = xx - rr*width;
    double x_r = xx + (1.0-rr)*width;
    if (step_out)
    {
        while (logdist(x_l) > log_uprime)
            x_l -= width;
        while (logdist(x_r) > log_uprime)
            x_r += width;
    }

    while (true)
    {
        xprime = dist(gen)*(x_r - x_l) + x_l;
        log_Px = logdist(xprime);
        if (log_Px > log_uprime)
            return xprime;
        else
        {
            if (xprime > xx)
                x_r = xprime;
            else if (xprime < xx)
                x_l = xprime;
            else
                throw 1;
        }
    }
}

#endif
