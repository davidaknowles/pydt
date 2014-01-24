#define _USE_MATH_DEFINES
#include <iostream>
#include <list>
#include <cstdlib>
#include <math.h>
//#include </home/dak33/boost_1_47_0/boost/lambda/lambda.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include "TreeNode.hpp"
#include <string>
#include <limits>
#include <algorithm>

using namespace std;

template <typename T>
T ListGet(list<T> l, int i)
{
    int counter=0;
    for (typename list<T>::iterator it=l.begin(); it!=l.end(); it++)
        if (counter==i)
            return *it;
        else
            counter++;
    throw 1;
}

// build a tree from the coalescent, return the root
TreeNode* KingmansCoalescent(int numLeaves, Settings &settings)
{
    list<TreeNode*> current;
    for (int i=0; i<numLeaves; i++)
        current.push_front(new TreeNode(0,true));
    for (int i=0; i<numLeaves-1; i++)
    {
        TreeNode *t = new TreeNode(-(i+1),false);
        int l = settings.iRand(0,numLeaves-i-1);
        TreeNode* l_node = ListGet<TreeNode*>(current,l);
        current.remove(l_node);
        int j = settings.iRand(0,numLeaves-i-2);
        TreeNode* j_node = ListGet<TreeNode*>(current,j);
        current.remove(j_node);
        t->children.push_back(l_node);
        t->children.push_back(j_node);
        current.push_back(t);
    }
    return *current.begin(); // should be the only element!
}
