
#ifndef VECTOR_EXTENSIONS_H
#define VECTOR_EXTENSIONS_H

#include <iostream>
#include <iomanip>
#include <boost/numeric/ublas/vector.hpp>
#include <list>

using namespace boost::numeric::ublas;

std::ostream &operator<<(std::ostream &output, const vector<double> &v)
{
    output << std::setprecision(10) ;
    for (int i=0; i<v.size(); i++)
        output << v[i] << " ";
    return output;
}

vector<double> &operator*=(vector<double> &left, const vector<double> &right)
{
    for (int i=0; i<left.size(); i++)
        left[i]*=right[i];
    return left;
}

vector<double> operator*(const vector<double> &left, const vector<double> &right)
{
    vector<double> result(left);
    result*=right;
    return result;
}

vector<double> &operator*=(vector<double> &left, double right)
{
    for (int i=0; i<left.size(); i++)
        left[i]*=right;
    return left;
}

vector<double> operator*(const vector<double> &left, double right)
{
    vector<double> result(left);
    result*=right;
    return result;
}

vector<double> &operator/=(vector<double> &left, const vector<double> &right)
{
    for (int i=0; i<left.size(); i++)
        left[i]/=right[i];
    return left;
}

vector<double> operator/(const vector<double> &left, const vector<double> &right)
{
    vector<double> result(left);
    result/=right;
    return result;
}

vector<double> &operator/=(vector<double> &left, double right)
{
    for (int i=0; i<left.size(); i++)
        left[i]/=right;
    return left;
}

vector<double> operator/(const vector<double> &left, double right)
{
    vector<double> result(left);
    result/=right;
    return result;
}

vector<double> &operator+=(vector<double> &left, const double &right)
{
    for (int i=0; i<left.size(); i++)
        left[i]+=right;
    return left;
}

vector<double> &operator+=(vector<double> &left, const vector<double> &right)
{
    for (int i=0; i<left.size(); i++)
        left[i]+=right[i];
    return left;
}

vector<double> operator+(const vector<double> &left, const double &right)
{
    vector<double> result(left);
    result+=right;
    return result;
}

vector<double> &operator-=(vector<double> &left, const double &right)
{
    for (int i=0; i<left.size(); i++)
        left[i]-=right;
    return left;
}

vector<double> operator-(const vector<double> &left, const double &right)
{
    vector<double> result(left);
    result-=right;
    return result;
}

vector<double> &operator&=(vector<double> &left, const double right)
{
    for (int i=0; i<left.size(); i++)
        left[i]=right;
    return left;
}

vector<double> applyFunc(const vector<double> &right, double (*func)(double))
{
    vector<double> result = right;
    for (int i=0; i<right.size(); i++)
        result[i]=func(right[i]);
    return result;
}

double sumVector(const vector<double> x)
{
    double result = 0.0;
    for (int i=0; i<x.size(); i++)
        result+=x[i];
    return result;
}

double min(const vector<double> x)
{
    double result = x[0];
    for (int i=1; i<x.size(); i++)
        if (x[i]<result)
            result=x[i];
    return result;
}

double max(const vector<double> x)
{
    double result = x[0];
    for (int i=1; i<x.size(); i++)
        if (x[i]>result)
            result=x[i];
    return result;
}

vector<double> listToVector(std::list<double> x)
{
    vector<double> result;
    result.resize(x.size());
    int j=0;
    for (std::list<double>::iterator i=x.begin(); i != x.end(); ++i)
    {
        result[j]=*i;
        j++;
    }
    return result;
}

double logSumExp(vector<double> x)
{
    double c = max(x);
    x -= c;
    vector<double> expon = applyFunc(x, exp);
    double s = sumVector(expon);
    return log(s)+c;
}

vector<double> softmax(vector<double> x)
{
    vector<double> result(x.size());
    double lse = logSumExp(x);
    for (int i=0; i<x.size(); i++)
    {
        result[i]=exp(x[i]-lse);
    }
    return result;
}

#endif
