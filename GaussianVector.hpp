#define _USE_MATH_DEFINES

#ifndef GAUSSIAN_VECTOR_H
#define GAUSSIAN_VECTOR_H

#include <math.h>
#include <boost/numeric/ublas/vector.hpp>
#include "VectorExtensions.hpp"
#include <numeric>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>

using namespace boost::numeric::ublas;

double reciprocal(double x)
{
    return 1.0 / x;
}

// Represents a diagonal multivariate Gaussian
class GaussianVector {

public:
    vector<double> MeanTimesPrecision;
    vector<double> Precision;
    bool isPoint; // whether this is a point mass, in which case MeanTimesPrecision stores the point
    static double lnSqrt2Pi;

    GaussianVector()
    {
    }

    GaussianVector(int D)
    {
        resize(D);
    }

    GaussianVector(vector<double> mp, vector<double> p)
    {
        MeanTimesPrecision = mp;
        Precision = p;
        isPoint = false;
    }

    int size()
    {
        return MeanTimesPrecision.size();
    }

    static vector<double> GaussianLogProb(vector<double> x, vector<double> mean, vector<double> variance)
    {
        vector<double> diff = x - mean;

        return -0.5 * (applyFunc(variance, log) + diff * diff / variance) - lnSqrt2Pi;
    }

    vector<double> GetLogProb(vector<double> x)
    {
        return GaussianLogProb(x, GetMean(), GetVariance());
    }

    double GetLogNormalizer()
    {
        if (isPoint)
            return 0.0;
        vector<double> terms = MeanTimesPrecision * MeanTimesPrecision / Precision - applyFunc(Precision, &log);
        return size() * lnSqrt2Pi + 0.5 * std::accumulate( terms.begin(), terms.end(), 0.0 );
    }

    void SetPoint(vector<double> point)
    {
        MeanTimesPrecision = point;
        isPoint = true;
    }


    vector<double> GetMean()
    {
        return isPoint ? MeanTimesPrecision : MeanTimesPrecision / Precision;
    }

    vector<double> GetVariance()
    {
        return isPoint ? MeanTimesPrecision * 0.0 : applyFunc(Precision, &reciprocal);
    }

    void SetMeanAndVariance(vector<double> mean, vector<double> variance)
    {
        isPoint = std::accumulate( variance.begin(), variance.end(), 0.0 ) == 0.0;
        if (isPoint)
            MeanTimesPrecision = mean;
        else
        {
            Precision = applyFunc(variance, &reciprocal);
            MeanTimesPrecision = mean * Precision;
        }
    }

    vector<double> Sample(boost::random::mt19937 rng)
    {
        boost::normal_distribution<double> norm;
        vector<double> result(size());
        for (int i=0; i<size(); i++)
            result[i] = norm(rng)/sqrt(Precision[i]) + MeanTimesPrecision[i]/Precision[i];
        return result;
    }

    static vector<double> Sample(boost::random::mt19937& rng, vector<double> mean, vector<double> variance)
    {
        boost::normal_distribution<double> norm;
        vector<double> result(mean.size());
        for (int i=0; i<mean.size(); i++)
            result[i] = norm(rng)*sqrt(variance[i]) + mean[i];
        return result;
    }

    void resize(int D)
    {
        MeanTimesPrecision.resize(D);
        Precision.resize(D);
    }

    void SetToUniform(int D)
    {
        MeanTimesPrecision.resize(D);
        MeanTimesPrecision &= 0.0;
        Precision.resize(D);
        Precision &= 0.0;
        isPoint = false;
    }

    GaussianVector &operator=(const GaussianVector& right)
    {
        MeanTimesPrecision = right.MeanTimesPrecision;
        Precision = right.Precision;
        isPoint = right.isPoint;
        return *this;
    }

    GaussianVector &operator*=(const GaussianVector& right)
    {
        if (isPoint)
        {
            return *this;
        }
        else if (right.isPoint)
        {
            *this = right;
        }
        else
        {
            MeanTimesPrecision += right.MeanTimesPrecision;
            Precision += right.Precision;
        }
        return *this;
    }

    GaussianVector &operator/=(const GaussianVector& right)
    {
        if (isPoint)
        {
            return *this;
        }
        else if (right.isPoint)
        {
            *this = right;
        }
        else
        {
            MeanTimesPrecision -= right.MeanTimesPrecision;
            Precision -= right.Precision;
        }
        return *this;
    }

    GaussianVector operator*(const GaussianVector& right)
    {
        GaussianVector result = *this;
        result *= right;
        return result;
    }

    GaussianVector operator/(const GaussianVector& right)
    {
        GaussianVector result = *this;
        result /= right;
        return result;
    }

    std::string str() __attribute__ ((noinline))
    {
        std::stringstream output;
        output << this;
        return output.str();
    }
};

double GaussianVector::lnSqrt2Pi = .5 * log( 2.0*M_PI );

std::ostream &operator<<(std::ostream &output, GaussianVector &v)
{
    output << "N(" << v.GetMean() << ", " << v.GetVariance() << ")";
    return output;
}

#endif
