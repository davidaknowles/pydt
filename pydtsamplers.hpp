#define _USE_MATH_DEFINES

#ifndef PYDTSAMPLERS_H
#define PYDTSAMPLERS_H

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <list>
#include <iterator>
#include <cstdlib>
#include <math.h>
#include <time.h>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/tokenizer.hpp>

#include <string>
#include <limits>
#include <algorithm>
#include <assert.h>

#include "slice_sample.hpp"
#include "Tree.hpp"

using namespace std;

// Functor used in slice sampling the hyperparameters
struct thetaHelper
{
    BaseSettings* settings;
    Tree* tree;

    thetaHelper( BaseSettings* _settings, Tree* _tree ) : settings(_settings), tree(_tree) {};

    double operator()(double theta)
    {
        if (theta < - 2.0 * settings->alpha)
            return log(0.0);
        settings->theta=theta;
        return tree->LogEvidenceStructure(*settings) - theta;
    }
};

// Functor used in slice sampling the hyperparameters
struct alphaHelper
{
    BaseSettings* settings;
    Tree* tree;

    alphaHelper( BaseSettings* _settings, Tree* _tree ) : settings(_settings), tree(_tree) {};

    double operator()(double alpha)
    {
        if ((settings->theta < - 2.0 * alpha) || (alpha < 0.0) || (alpha >= 1.0))
            return log(0.0);
        settings->alpha=alpha;
        return tree->LogEvidenceStructure(*settings);
    }
};

// Functor used in slice sampling the hyperparameters
struct cHelper
{
    BaseSettings* settings;
    Tree* tree;

    cHelper( BaseSettings* _settings, Tree* _tree ) : settings(_settings), tree(_tree) {};

    double operator()(double c)
    {
        if (c < 0.0)
            return log(0.0);
        settings->c=c;
        return tree->LogEvidenceStructure(*settings) - c;
    }
};

void SampleHypers(Tree* tree, BaseSettings& settings)
{
    if (settings.sampleTheta)
    {
        thetaHelper th(&settings,tree);
        settings.theta = slice_sample(settings.theta, th, settings.gen);
    }

    if (settings.sampleAlpha)
    {
        alphaHelper ah(&settings,tree);
        settings.alpha = slice_sample(settings.alpha, ah, settings.gen);
    }

    if (settings.sampleC)
    {
        cHelper ch(&settings,tree);
        settings.c = slice_sample(settings.c, ch, settings.gen);
    }
}

/* The simple Metropolis Hastings sampler. 
 Inputs: 
    - tree: initial tree
    - x_star: test data
    - settings: algorithm settings and hyperparameters
    - iterations: number of MH steps to take
    - filename: log file name
    - newickfile: file to store learnt structure in newick format
 Output: new tree structure */ 
Tree* MH(Tree* tree, list<boost::numeric::ublas::vector<double> > &x_star, BaseSettings& settings, int iterations, const char* filename,  const char* newickfile = "")
{
    remove(filename);

    long accepts = 0 ;
    double ml = tree->MarginalLikelihood(settings);
    long subtreeSizeCounter = 0 ;
    clock_t start_time = clock();
    for (int i=0; i<iterations; i++)
    {
        tree->countLeaves(true);
        Tree* newTree = tree->DeepCopy();
        TreeNode* subtree = newTree->GetUniformRandomSubtree(settings);

        double oldAttachmentProb = newTree->GetProbOfAttachSubtree(settings, subtree);

        subtreeSizeCounter += subtree->countLeaves();
        newTree->DetachSubtree(subtree);

        bool success=newTree->AttachSubtree(settings,subtree);
        if (success)
        {
            double newAttachmentProb = newTree->GetProbOfAttachSubtree(settings,subtree);
            assert(tree->root->numDescendants==newTree->root->numDescendants);
            double newMl = newTree->MarginalLikelihood(settings);
            double a = newMl-ml+oldAttachmentProb-newAttachmentProb;
            if (exp(a) > settings.fRand())
            {
                delete tree;
                tree=newTree;
                ml=newMl;
                accepts++;
            }
            else
                delete newTree;

        }
        else
            delete newTree;

        SampleHypers(tree, settings);

        int printEvery = 50;
        if (i % printEvery == 1)
        {
            double logPred=0.0;//, logPred2=0.0;
            if (!x_star.empty())
            {
                boost::numeric::ublas::vector<double> pred = tree->Prediction(x_star, settings);
                logPred=sumVector(applyFunc(pred,log));
            }

            int numInternalNodes = tree->numInternalNodes();
            double averageBranchingFactor = tree->averageBranchingFactor(); 
            int maxBranchingFactor = tree->root->maxBranchingFactor();

            cout << "It " << i << " time(s) " << (double)(clock()-start_time)/(double)CLOCKS_PER_SEC << " ML: " << ml << " pred " << logPred <<  " theta: " << settings.theta << " alpha: " << settings.alpha << " c: " << settings.c << " internal nodes " << numInternalNodes << " brachingFactor " << averageBranchingFactor << "(max " << maxBranchingFactor << ") ar: " << (double)accepts/(double)printEvery << endl;
            accepts=0;
            ofstream outputfile;

            if (newickfile != "")
            {
                outputfile.open(newickfile,ios::out);
                outputfile  << tree->root->newick(0.0) << endl;
                outputfile.close();
            }

            outputfile.open(filename,ios::app);
            outputfile  << i << " " << (double)(clock()-start_time)/(double)CLOCKS_PER_SEC << " " << ml  << " " << logPred << " " << numInternalNodes <<  " " << averageBranchingFactor << " " << maxBranchingFactor << " " << settings.theta << " " << settings.alpha << " " << settings.c << endl;
            outputfile.close();
        }
    }
    return tree;
}

/* The slice sampler. 
 Inputs: 
    - tree: initial tree
    - x_star: test data
    - settings: algorithm settings and hyperparameters
    - iterations: number of slice sampling moves to make
    - filename: log file name
    - newickfile: file to store learnt structure in newick format
 Output: new tree structure */ 
Tree* SliceSample(Tree* tree, list<boost::numeric::ublas::vector<double> > &x_star, BaseSettings& settings, int iterations, const char* filename, const char* newickfile="")
{
    int numLeaves = tree->countLeaves();
    remove(filename);
    double ml = tree->MarginalLikelihood(settings);
    long subtreeSizeCounter = 0 ;
    long attemptsCounter = 0;
    clock_t start_time = clock();
    for (int i=0; i<iterations; i++)
    {
        // slice height in log space
        double logu = ml + log(settings.fRand());

        // events store rejected points
        tree->ClearEvents();

        // choose a subtree uniformly
        TreeNode* subtree = tree->GetUniformRandomSubtree(settings);

        // keep track of the size of detached subtrees
        subtreeSizeCounter += subtree->countLeaves();

        // Store where this subtree was originally attached: need this to shrink the slice
        BranchPoint originalPosition;
        tree->DetachSubtree(subtree,&originalPosition);

        if (originalPosition.isExistingBranchPoint)
        {
            logu -= log(settings.nodeWidth); // we represent atoms of mass m as rectangles of width nodeWidth and height m/nodeWidth
            originalPosition.time = settings.fRand(0.0, settings.nodeWidth); // define a pretend position in the node "slice"
        }

        // Defines the upper limit of the slice (note that we start with the entire tree
        BranchPoint topOfSlice(tree->zero,tree->root,0.0,false);
        bool firstTime=true;
        int attempts = 0;
        double oldSlice = 1e10; // used to check the slice always gets smaller
        while (true) // shrink
        {
            // if not the first iteration then the subtree will have been attached and needs removing
            if (!firstTime)
                tree->DetachSubtree(subtree);
            firstTime=false;
            tree->countLeaves(true);
            double length=0.0;
            int numNodes=0;

            // this does the shrinking as well
            BranchPoint bp = tree->SampleFromSlice(settings,topOfSlice,originalPosition,subtree->time,length,numNodes);

            if (!settings.multifurcating())
                assert( !bp.isExistingBranchPoint ); 

            assert( length <= oldSlice );
            oldSlice = length;
            assert( tree->zero->children.size() == 1);

            bp.AttachSubtree(subtree, settings);

            tree->root=tree->zero->children.front();
            assert( tree->zero->children.size() == 1);

            assert( numLeaves == tree->countLeaves() );

            double newMl = tree->MarginalLikelihood(settings);
            assert(!isnan(newMl));

            attempts++;
            assert(attempts < 1000);

            if (bp.time == originalPosition.time || (newMl - (bp.isExistingBranchPoint ? log(settings.nodeWidth) : 0.0) >= logu))
            {
                ml=newMl;
                break;
            }
        }
        attemptsCounter += attempts;

        SampleHypers(tree, settings);

        int printEvery = 10;
        if (i % printEvery == 1)
        {

            double logPred=0.0;//, logPred2=0.0;
            if (!x_star.empty())
            {
                boost::numeric::ublas::vector<double> pred = tree->Prediction(x_star, settings);
                logPred=sumVector(applyFunc(pred,log));
            }

            int numInternalNodes = tree->numInternalNodes();
            double averageBranchingFactor = tree->averageBranchingFactor(); 
            int maxBranchingFactor = tree->root->maxBranchingFactor();

            cout << "It " << i << " time(s) " << (double)(clock()-start_time)/(double)CLOCKS_PER_SEC << " ML: " << ml << " pred " << logPred <<  " theta: " << settings.theta << " alpha: " << settings.alpha << " c: " << settings.c << " # internal " << numInternalNodes << " branchingFactor " << averageBranchingFactor << "(max " << maxBranchingFactor << ") attempts " << ((double)attemptsCounter/(double)printEvery)<< endl;
            attemptsCounter = 0;

            ofstream outputfile;
            if (newickfile != "")
            {
                outputfile.open(newickfile,ios::out);
                outputfile  << tree->root->newick(0.0) << endl;
                outputfile.close();
            }

            outputfile.open(filename,ios::app);
            outputfile  << i << " " << (double)(clock()-start_time)/(double)CLOCKS_PER_SEC << " " << ml  << " " << logPred << " " << numInternalNodes <<  " " << averageBranchingFactor << " " << maxBranchingFactor << " " <<settings.theta << " " << settings.alpha << " " << settings.c << endl;
            outputfile.close();
        }


    }
    return tree;
}

// -------------- EXPERIMENTAL -----------------------
Tree* PoissonProcessSampler(Tree* tree, BaseSettings& settings, int iterations)
{
    tree->InstantiateLocations(settings);

    for (int i=0; i<iterations; i++)
    {
        cout << "Current tree: " << tree->newickWithTime() << endl;
        TreeNode* subtree = tree->GetRandomLeafAndSampleRST(settings); // subtree can only be a leaf for now

        tree->DetachSubtree(subtree);
        tree->EvaluatePotentialBranchPoints(subtree, settings);
        // note: if we calculate the ML here it will overwrite all the location samples
    }
    return tree;
}

// -------------- EXPERIMENTAL -----------------------
Tree* PoissonProcessSamplerWrapper(Tree* tree, BaseSettings& settings, int outerIts, int innerIts, const char* filename,  list<boost::numeric::ublas::vector<double> > &test)
{

    for (int i=0; i<outerIts; i++)
    {
        tree=PoissonProcessSampler(tree, settings, innerIts);
        double logPredictive=0;
        //double logPredictive = sumVector(tree->PredictionSimple(test, settings));
        double ml=tree->MarginalLikelihood(settings);
//     cout << "Outer it " << i << " ML: " << ml << " pred: " << logPredictive << endl;
        ofstream outputfile;
        outputfile.open(filename,ios::app);
        outputfile << tree->root->time << endl;
        outputfile.close();
    }
    //tree->OutputTree("tree.txt");
    return tree;
}


Tree* BuildTree(list<TreeNode*> &leaves, BaseSettings& settings)
{
    TreeNode* root=new TreeNode(0.5,false);
    list<TreeNode*>::iterator i=leaves.begin();
    root->children.push_back((*i)->DeepCopy());
    i++;
    root->children.push_back((*i)->DeepCopy());
    Tree *tree = new Tree(root,settings);
    for (i++; i != leaves.end(); i++)
    {
        tree->InstantiateLocations(settings);
        tree->SampleTPrior(settings);
        tree->EvaluatePotentialBranchPoints((*i)->DeepCopy(), settings);
    }
    return tree;
}

#endif
