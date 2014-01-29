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

	double mlRemaining=newTree->MarginalLikelihoodOnly(settings); // get up and down messages
	BranchPoint where; 
        bool success=newTree->AttachSubtree(settings,subtree,where);
        if (success)
        {
	  double newMlCheck; 
	  GaussianVector temp(settings.D);
	  double subtreeML = subtree->bpSweepUp2(temp,false,settings.D); 
	  if (where.isExistingBranchPoint){
	    subtree->Msg_Normal_ParentLocation = SampleAverageConditional(temp,subtree->time - where.parent->time); 
	    double diff= where.parent->local_factor(false,settings.D) - where.parent->local_factor(true,settings.D); 
	    newMlCheck = diff + subtreeML + mlRemaining; 
	    cout << " newMlcheck: " << newMlCheck << " newMl " << newTree->MarginalLikelihoodOnly(settings) << endl;
	  } else {
	    // remove contribution from p_i in pi_l, imagining it was root
	    // add contribution from p_i in newtree, with p_l as root
	    // add contribution from p_l in newtree, as root
	    TreeNode* new_node= where.parent->children.back(); 
	    new_node->Msg_Normal_ParentLocation = where.child->Msg_Normal_ParentLocation; 
	    double old_pi_contribution; 
	    double new_pi_contribution; 
	    if (where.parent->time==0.0){ // attached directly below the root! 
	      boost::numeric::ublas::vector<double> zeros(settings.D);
	      zeros &= 0.0;
	      old_pi_contribution = sumVector(where.child->Msg_Normal_ParentLocation.GetLogProb(zeros)); 
	      new_pi_contribution = 0.0; 
	    } else {
	      old_pi_contribution = where.parent->local_factor(false,settings.D); 
	      new_pi_contribution = where.parent->local_factor(true,settings.D); 
	    }
	    GaussianVector wp_to_norm = where.parent->Marg_Location / where.child->Msg_Normal_ParentLocation; 
	    new_node->Msg_Normal_Location = SampleAverageConditional(wp_to_norm, new_node->time - where.parent->time) ;
	    GaussianVector wc_to_norm = where.child->Marg_Location / where.child->Msg_Normal_Location; 
	    where.child->Msg_Normal_ParentLocation = SampleAverageConditional(wc_to_norm, where.child->time - new_node->time); 
	    subtree->Msg_Normal_ParentLocation = SampleAverageConditional(temp,subtree->time - new_node->time); 
	    double pl_contribution = new_node->local_factor(false,settings.D); 
	    double myML = mlRemaining + subtreeML + new_pi_contribution + pl_contribution - old_pi_contribution; 
	    if (isnan(myML)) throw 1; 
	    assert( abs(myML - newTree->MarginalLikelihoodOnly(settings)) < 0.001); 
	    cout << " myML " << myML << " true " << newTree->MarginalLikelihoodOnly(settings) << endl;
	  }

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
    long subtreeSizeCounter = 0 ;
    double current_ml=tree->MarginalLikelihoodOnly(settings);
    double current_struct_prior=tree->LogEvidenceStructureOnly(settings); 
    double current_time_prior=tree->LogEvidenceTimes(settings); 
    tree->SetupParents(); 
    long attemptsCounter = 0;
    clock_t start_time = clock();
    
    // TODO: should NEVER have to do full BP sweep. 
    // When a subtree is detached you already have it's up messages, and only the messages _up_ to the root change from the detachment position changes
    // When the subtree is reattached, the messages _down_ the subtree change, and _up_ to the root

   

    for (int i=0; i<iterations; i++)
    {
      // tree->root->CheckConsistency();

      double current_log_joint = current_ml + current_time_prior + current_struct_prior; 

      // slice height in log space
      double logu = current_log_joint + log(settings.fRand());
      
      // events store rejected points
      tree->ClearEvents();
      
      // choose a subtree uniformly
      TreeNode* subtree = tree->GetUniformRandomSubtree(settings);

        // keep track of the size of detached subtrees
        subtreeSizeCounter += subtree->countLeaves();

	tree->SetupParents(true); 
        // Store where this subtree was originally attached: need this to shrink the slice
	double ml_change,time_prior_change,struct_prior_change; 
	BranchPoint originalPosition = subtree->Detach(settings,ml_change,time_prior_change,struct_prior_change); 
	tree->root=tree->zero->children.front();

	tree->root=tree->zero->children.front();
	current_ml += ml_change;
	current_time_prior += time_prior_change; 
	current_struct_prior += struct_prior_change; 
	// upwards sweep on the subtree
	GaussianVector msg_from_subtree = subtree->Marg_Location / subtree->Msg_Normal_Location; 
	
	if (settings.debug){
	  double mlRemaining = tree->MarginalLikelihoodOnly(settings); 
	  GaussianVector msg_from_subtree_check(settings.D); 
	  double subtreeML = subtree->bpSweepUp2(msg_from_subtree_check,false,settings.D); 
	  double timesSubtree = subtree->LogEvidenceTimes(0.0,settings) - subtree->LocalEvidenceTimes(0.0,settings); 
	  double timesRemaining = tree->LogEvidenceTimes(settings); 
	  double structureSubtree = subtree->LogEvidenceStructureOnly(settings); 
	  double structureRemaining = tree->LogEvidenceStructureOnly(settings); 

	  assert( abs( current_ml - (mlRemaining + subtreeML) ) < 0.001); 
	  assert( abs( current_struct_prior - (structureSubtree + structureRemaining) ) < 0.001); 
	  assert( abs( current_time_prior - (timesSubtree + timesRemaining) ) < 0.001);
	  assert( abs( msg_from_subtree.GetMean()[0] - msg_from_subtree_check.GetMean()[0] ) < 0.001); 
	}

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
            firstTime=false;
	    if (settings.debug)
	      tree->countLeaves(true);
            double length=0.0;
            int numNodes=0;

            // this does the shrinking as well
            BranchPoint bp = tree->SampleFromSlice(settings,topOfSlice,originalPosition,subtree->time,length,numNodes);
	    if (bp.isExistingBranchPoint)
	      assert( bp.child->time > 0.0); 
            if (!settings.multifurcating())
	      assert( !bp.isExistingBranchPoint ); 

            assert( length <= oldSlice );
            oldSlice = length;
            assert( tree->zero->children.size() == 1);
	    // fast calculation [note we want to do this BEFORE attaching]
	    double new_ml = current_ml; 
	    double new_time_prior = current_time_prior; 
	    double new_struct_prior = current_struct_prior; 

	    if (bp.isExistingBranchPoint){
	      TreeNode* atnode=bp.child; 
	      double old_contrib = atnode->local_factor(false,settings.D);
	      double old_prior_contrib = atnode->LogEvidenceStructureOnlyUp(settings,0);
	      double old_times_contrib = atnode->LogEvidenceTimesOnlyUp(settings,0); 
	      atnode->children.push_back(subtree); 
	      subtree->Msg_Normal_ParentLocation = SampleAverageConditional(msg_from_subtree,subtree->time - atnode->time); 
	      double new_contrib = atnode->local_factor(false,settings.D);
	      new_ml += new_contrib - old_contrib; 
	      new_struct_prior += atnode->LogEvidenceStructureOnlyUp(settings,subtree->numDescendants) - old_prior_contrib; 
	      new_time_prior += atnode->LogEvidenceTimesOnlyUp(settings,subtree->numDescendants) - old_times_contrib; 
	      new_time_prior += subtree->LocalEvidenceTimes(atnode->time,settings); // - subtree->LocalEvidenceTimes(originalPosition.time,settings); 
	      atnode->children.pop_back();
	    } else {
	      double old_pi_contribution,new_pi_contribution;
	      if (bp.parent->time==0.0){
		boost::numeric::ublas::vector<double> zeros(settings.D);
		zeros &= 0.0;
		old_pi_contribution = sumVector(bp.child->Msg_Normal_ParentLocation.GetLogProb(zeros)); 
		new_pi_contribution = 0.0; 
	      } else {
		old_pi_contribution = bp.parent->local_factor(false,settings.D); 
		// want parent WITHOUT bp.child included
		new_pi_contribution = bp.parent->local_factor(false,settings.D,bp.child);
	      }
	      TreeNode* new_node=new TreeNode(bp.time,false); 
	      GaussianVector wp_to_norm = bp.parent->Marg_Location / bp.child->Msg_Normal_ParentLocation; 
	      new_node->Msg_Normal_Location = SampleAverageConditional(wp_to_norm, new_node->time - bp.parent->time) ;
	      GaussianVector wc_to_norm = bp.child->Marg_Location / bp.child->Msg_Normal_Location; 
	      new_node->children.push_back(subtree); 
	      subtree->Msg_Normal_ParentLocation = SampleAverageConditional(msg_from_subtree,subtree->time - new_node->time);
	      TreeNode* bp_child_copy = new TreeNode(bp.child->time, bp.child->isLeaf); 
	      bp_child_copy->numDescendants = bp.child->numDescendants;
	      bp_child_copy->Msg_Normal_ParentLocation = SampleAverageConditional(wc_to_norm, bp.child->time - new_node->time); 
	      new_node->children.push_back(bp_child_copy); 
	      double pl_contribution = new_node->local_factor(false,settings.D); 
	      new_ml += new_pi_contribution + pl_contribution - old_pi_contribution; 
	      new_node->numDescendants = subtree->numDescendants + bp.child->numDescendants;
	      new_time_prior += new_node->LocalEvidenceTimes(bp.parent->time,settings); 
	      new_struct_prior += new_node->LocalEvidenceStructure(settings); 
	      new_time_prior += subtree->LocalEvidenceTimes(new_node->time,settings); // - subtree->LocalEvidenceTimes(originalPosition.time,settings); 
	      new_time_prior += bp.child->LocalEvidenceTimes(new_node->time,settings) - bp.child->LocalEvidenceTimes(bp.parent->time,settings); 
	      new_struct_prior -= bp.parent->LogEvidenceStructureOnlyUp(settings,0);
	      new_time_prior -= bp.parent->LogEvidenceTimesOnlyUp(settings,0); 
	      bp.child->numDescendants += subtree->numDescendants; 
	      new_struct_prior += bp.parent->LogEvidenceStructureOnlyUp(settings,subtree->numDescendants); 
	      new_time_prior += bp.parent->LogEvidenceTimesOnlyUp(settings,subtree->numDescendants); 
	      bp.child->numDescendants -= subtree->numDescendants; 
	      delete bp_child_copy; 
	      delete new_node; 
	    }
	    
	    assert(!isnan(new_ml));

	    double new_log_joint = new_ml + new_struct_prior + new_time_prior; 
	    
            attempts++;
            assert(attempts < 1000);

            if (bp.time == originalPosition.time || (new_log_joint - (bp.isExistingBranchPoint ? log(settings.nodeWidth) : 0.0) >= logu))
            {
	      if (settings.debug){
		tree->countLeaves(true);
		subtree->countLeaves(true,false); 
	      }

	      bp.AttachSubtree(subtree, settings);
	      tree->root=tree->zero->children.front();
	      
	      if (settings.debug){
		assert(numLeaves == tree->countLeaves(true)); 
		tree->SetupParents(true);  // check parents are correctly setup
		assert( abs( new_ml - tree->MarginalLikelihoodOnly(settings)) < 0.001); 
		assert( abs(new_struct_prior - tree->LogEvidenceStructureOnly(settings)) < 0.001); 
		assert( abs(new_time_prior - tree->LogEvidenceTimes(settings)) < 0.001); 
	      }
	      current_time_prior = new_time_prior; 
	      current_struct_prior = new_struct_prior; 
	      current_ml = new_ml ; 

	      assert( tree->zero->children.size() == 1);
	      break;
            }
        }
        attemptsCounter += attempts;

	int hypers_every = 50; 
	if (i % hypers_every == (hypers_every-1)){
	  SampleHypers(tree, settings);
	  settings.HypersChanged(); 
	  current_struct_prior=tree->LogEvidenceStructureOnly(settings); 
	  current_time_prior=tree->LogEvidenceTimes(settings); 
	  current_ml = tree->MarginalLikelihoodOnly(settings); // just to refresh messages
	}
	  

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
	    
            cout << "It " << i << " time(s) " << (double)(clock()-start_time)/(double)CLOCKS_PER_SEC << " ML: " << current_ml + current_time_prior + current_struct_prior << " pred " << logPred <<  " theta: " << settings.theta << " alpha: " << settings.alpha << " c: " << settings.c << " # internal " << numInternalNodes << " branchingFactor " << averageBranchingFactor << "(max " << maxBranchingFactor << ") attempts " << ((double)attemptsCounter/(double)printEvery)<< endl;
            attemptsCounter = 0;

            ofstream outputfile;
            if (newickfile != "")
            {
                outputfile.open(newickfile,ios::out);
                outputfile  << tree->root->newick(0.0) << endl;
                outputfile.close();
            }

            outputfile.open(filename,ios::app);
            outputfile  << i << " " << (double)(clock()-start_time)/(double)CLOCKS_PER_SEC << " " << current_ml + current_time_prior + current_struct_prior << " " << logPred << " " << numInternalNodes <<  " " << averageBranchingFactor << " " << maxBranchingFactor << " " <<settings.theta << " " << settings.alpha << " " << settings.c << endl;
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
