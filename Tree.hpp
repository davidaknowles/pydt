#define _USE_MATH_DEFINES

#ifndef TREE_H
#define TREE_H

#include <iostream>
#include <list>
#include <map>
#include <cstdlib>
#include <math.h>
#include <string>
#include <limits>
#include <numeric>
#include <assert.h>
#include "TreeNode.hpp"

using namespace std;

// Represents a tree structure.
class Tree {

public:
    // The root of the tree
    TreeNode* root;

    // The "zero" node which lives at (0,0), useful as the parent of root.
    TreeNode* zero;

    // Make a deep copy of the tree by recursing down the structure.
    // Note that messages are not copied, but data at the leaves is copied
    Tree* DeepCopy()
    {
        Tree* result = new Tree(); // the new tree
        result->zero=zero->DeepCopy(); // recurse
        result->zero->Marg_Location=zero->Marg_Location; // Keep zero at zero!
        result->root=result->zero->children.front(); // set the root of the new root. Zero only has one child
        return result;
    }

    // Constructor
    Tree()
    {

    }

    // Better constructor: r is the root, settings
    Tree(TreeNode* r, BaseSettings& settings)
    {
        root=r;
        // Setup the zero node
        boost::numeric::ublas::vector<double> zeros(settings.D);
        zeros &= 0.0;

        zero = new TreeNode(0,false);
        zero->Marg_Location.SetPoint(zeros);
        zero->numDescendants=root->numDescendants;
        zero->children.push_back(root);
    }

    // Destructor. Recursively dealloc all memory for this tree
    ~Tree()
    {
        TreeNode::deleteSubtree(zero);
    }

    // Generate a tree copying the list of leaves provided, where all the leaves
    // attach to the root directly
    static Tree* FlatTree(list<TreeNode*> &leaves, BaseSettings& settings)
    {
        TreeNode* root = new TreeNode(0.5, false);
        root->numDescendants=leaves.size();
        //cout << "Num leaves: " << leaves.size() << endl;
        for (list<TreeNode*>::iterator i=leaves.begin(); i != leaves.end(); ++i)
        {
            root->children.push_back((*i)->DeepCopy());
        }
        return new Tree(root, settings);
    }

    static Tree* RandomTree(list<TreeNode*> &leaves, BaseSettings& settings)
    {
        // initially the root will be the first leaf!
        list<TreeNode*>::iterator i=leaves.begin();
        TreeNode* root = (*i)->DeepCopy();
        root->numDescendants=1;
        Tree *tree = new Tree(root,settings);
        for (i++; i != leaves.end(); i++)
        {
            tree->AddChild(settings, (*i)->DeepCopy());
            //cout << tree->newick() << endl;
        }
        //cout << "Num leaves: " << tree->countLeaves() << endl;
        return tree;
    }

    // Detach a subtree (defined by its root) from the tree
    void DetachSubtree(TreeNode* subtree, BranchPoint* originalPosition = NULL)
    {
        // This may change the root of the tree
        root=root->DetachSubtree(subtree,zero, originalPosition);
        // don't need to delete the old root, this is handled by the above recursion
        zero->children.clear(); // zero's child should be the new root
        zero->children.push_back(root);
        // the cached number of descendants for zero must be updated
        zero->numDescendants=root->numDescendants;
    }

    // Change the current root and delete the old one.
    void SetNewRoot(TreeNode *newRoot)
    {
        if (root==newRoot) // check these are not equal
            throw 1; // could ignore this but think I'd rather know about it!
        //cout << root->newick() << endl;
        delete root; // delete the old root
        root=newRoot;
        zero->children.clear(); // zero's child should be the new root
        zero->children.push_back(root);
    }

    // Allocate memory for messages prior to running BP
    void initialise(BaseSettings &settings)
    {
        zero->initialise(settings); // allocated messages recursively
        // Set zero vector
        boost::numeric::ublas::vector<double> zeros(settings.D);
        zeros &= 0.0;
        zero->Marg_Location.SetPoint(zeros);
    };

    // Add a new leaf by running the generative process
    // Returns what depth this leaf was attached at
    int AddChild(BaseSettings& s, TreeNode* leaf) {
        zero->numDescendants++;
        double At = exp(s.logDivergenceRateFactor(root->numDescendants))*(-log(s.fRand(0,1)));
        double td= s.invA(At);
        if (td>root->time) // don't diverge off this branch
            return root->AddChild(s,leaf);
        else // diverge
        {
            // create new internal node
            TreeNode * temp = new TreeNode(td,false);
            temp->children.push_back(root);
            temp->children.push_back(leaf);
            temp->numDescendants = root->numDescendants+1;
            root=temp;
            zero->children.clear(); // zero's child should be the new root
            zero->children.push_back(root);
            return 1;
        }
    };

    // Try to attach the subtree using the generative process
    // Return whether this was successful
  bool AttachSubtree(BaseSettings& s, TreeNode* subtree, BranchPoint &where, int maxTries = 100) {
        bool success=false; // whether we have been successful yet
        // try up to maxTries times to attach the subtree
        for (int j=0; j<maxTries && !success; j++)
        {
	  success = zero->AttachSubtree(s,subtree,where,true);
        }
        // The root may have changed, update this
        if (success)
            root=zero->children.front();
        return success;
    };

    double GetProbOfAttachSubtree(BaseSettings& s, TreeNode* subtree)
    {
        double logProb = (0.0-s.A(root->time)) * exp(s.singleTerm(root->numDescendants-subtree->numDescendants));
        bool success= root->GetProbOfAttachSubtree(s, subtree, logProb);
        if (!success)
            throw 1;
        return logProb;
    }

    // Check consistency of messages
    void CheckConsistency()
    {
        root->CheckConsistency();
    }

    void ClearEvents()
    {
        zero->ClearEvents();
    }


    // Output tree in newick format including divergence times
    string newickWithTime()
    {
        return root->newick(0.0);
    };

    //  Output tree in newick format including divergence times AND locations
    string newickTimeAndLocation()
    {
        return root->newick2(0.0);
    };

    // Output tree in newick format (structure only)
    string newick()
    {
        return root->newick();
    };

    // Count the number of leaves.
    // check: whether to throw an exception if numDescendants is incorrect
    int countLeaves(bool check = false, bool countAll = false)
    {
        return root->countLeaves(check, countAll);
    };

    int numInternalNodes()
    {
        return root->countLeaves(false, true) - root->countLeaves(false, false);
    };

    double averageBranchingFactor()
    {
        return (double)root->sumBranchingFactor() / (double)numInternalNodes(); 
    }

    // Sample data at the leaves
    void sampleData(BaseSettings &s)
    {
        boost::numeric::ublas::vector<double> zeros(s.D);
        zeros &= 0.0;
        zero->sampleData(zeros, 0.0, s);
    };

    // log evidence contribution of the locations, calculated EP style (probably can do this more straightforwardly!)
    // I've compared this to Infer.NET to make sure it's correct
    double LogEvidenceSimple(BaseSettings &s)
    {
      GaussianVector temp;         
      boost::numeric::ublas::vector<double> zeros(s.D);
      zeros &= 0.0;
      double ml = root->bpSweepUp2(temp,true,s.D); 
      return ml; // + sumVector(temp.GetLogProb(zeros)); 
    }

    // Calculate the marginal likelihood of this tree structure
    double MarginalLikelihood(BaseSettings &settings)
    {
        initialise(settings); // allocate memory for messages
        root->bpSweepUp(*zero); // belief propagation
        root->bpSweepDown(*zero);
	double logEvidence=root->LogEvidence(*zero); 
	double leCheck=LogEvidenceSimple(settings); 
	if (abs(logEvidence-leCheck) > 0.001){
	  cout << "old: " << logEvidence << " new: " << leCheck << endl; 
	  throw 1; 
	}
        return logEvidence + LogEvidenceStructure(settings);
    }

    double LogEvidenceStructure(BaseSettings &settings) // include P(times|structure)
    {
        return root->LogEvidenceStructure(*zero, settings) ;
    }

    double LogEvidenceStructureOnly(BaseSettings &settings) // P(structure), no times
    {
        return root->LogEvidenceStructureOnly(*zero, settings) ;
    }

    double LogEvidenceTimes(BaseSettings &settings)
    {
        return root->LogEvidenceTimes(*zero, settings) ;
    }

    // Calculate the marginal likelihood of this tree structure
    double MarginalLikelihoodOnly(BaseSettings &settings)
    {
        initialise(settings); // allocate memory for messages
        root->bpSweepUp(*zero); // belief propagation
        root->bpSweepDown(*zero);
        return root->LogEvidence(*zero);
    }

    // Calculate the marginal likelihood of this tree structure
    void InstantiateLocations(BaseSettings &settings)
    {
        initialise(settings); // allocate memory for messages
        root->bpSweepUp(*zero); // belief propagation
        root->bpSweepDown(*zero);
        root->sampleSweepUp(*zero,settings);
        //CheckConsistency();
    }

    // Return a random subtree of the tree, but NOT the root. Using settings.subtreeInflationFactor
    // you can control how likely larger subtrees are to be detached
    TreeNode* GetRandomSubtree(BaseSettings &settings)
    {
        return root->GetRandomSubtree(settings);
    }

    // Return a random leaf of the tree.
    // you can control how likely larger subtrees are to be detached
    TreeNode* GetRandomLeaf(BaseSettings &settings)
    {
        list<TreeNode*> leaves;
        root->GetAllLeaves(leaves);
        int index=settings.iRand(0,leaves.size()-1);
        return listIndexOf(leaves,index);
    }

    // Return a uniform at random subtree of the tree, but NOT the root
    TreeNode* GetUniformRandomSubtree(BaseSettings &settings)
    {
        list<TreeNode*> allNodes;
        root->GetAllNodes(allNodes);
        int index=settings.iRand(0,allNodes.size()-1);
        return listIndexOf(allNodes,index);
    }

    /* This function performs the main slice sampling logic. It builds a list of possible reattachment positions, represented by BranchPoint
    instances, and their associated lengths. This includes both positions on branches (then branchpoint.isExistingBranchPoint=false) or 
    at existing nodes (then branchpoint.isExistingBranchPoint=true). 
    Inputs: 
        - topOfSlice: the upper (earliest, smallest time) limit of the slice, 
        - originalPosition: where the subtree was originally attached, 
        - subtreeTime: time for the root of the deattached subtree
    Outputs: 
        - length: total length of the slice
        - numNodes: number of branchPoints in the slice (branches and nodes) 
        - return value: proposed attachment position (to be accepted or rejected) 
    
    Notes: The logic here gets a little complex. topOfSlice.time represents the upper limit of the slice. If topOfSlice is a branch
    (i.e. topOfSlice.isExistingBranchPoint=false) then topOfSlice.time is simply the time of the upper extent of the slice. 
     If topOfSlice is a node (i.e. topOfSlice.isExistingBranchPoint=true) then time is how far into the interval representing
    the top node the upper limit of the slice is. For branches, the parent node has a member "events", a sorted list of rejected
    attachment positions on that branch, so that the first element gives a lower limit of the slice. For nodes, the member "node_event" 
    records whether there has been a rejected reattachment at this node, in which case node_event > 0, with the value representing how 
    much of the interval corresponding to this node was used up (the initial interval being settings.nodeWidth).  
    */  
    BranchPoint SampleFromSlice(BaseSettings &settings, BranchPoint& topOfSlice, BranchPoint& originalPosition, double subtreeTime, double &length, int &numNodes, bool debug=false)
    {
        list<BranchPoint> nodes; // possible reattachment positions
        list<double> lengths; // lengths of these positions
        topOfSlice.onPath=true; // "onPath" is whether this node is on the path from the root to the original attachment position
        nodes.push_back(topOfSlice);

        // ->events contains rejected attachment positions
        if (topOfSlice.child->events.empty() || topOfSlice.isExistingBranchPoint) // then we should carry on down the tree
        {
            lengths.push_back(topOfSlice.isExistingBranchPoint ?
                              (settings.nodeWidth - topOfSlice.time - topOfSlice.child->node_event) :
                              (min(topOfSlice.child->time,subtreeTime)-topOfSlice.time));

            assert( lengths.back() > 0.0 );

            // should we add a potential attachment position at topOfSlice.child
            if ((topOfSlice.child->time < subtreeTime) && settings.multifurcating() && !topOfSlice.isExistingBranchPoint)
            {
                BranchPoint atnode=topOfSlice;
                atnode.isExistingBranchPoint=true;
                nodes.push_back(atnode);
                lengths.push_back(settings.nodeWidth - topOfSlice.child->node_event); // if topOfSlice is a node then topOfSlice.time is the amount of slice it has left
            }
            if (topOfSlice.child->node_event==0.0)
            {
                for (list<TreeNode*>::iterator i=topOfSlice.child->children.begin(); i != topOfSlice.child->children.end(); ++i)
                {
                    (*i)->GetBranchesInSlice(settings,nodes,lengths,*topOfSlice.child,originalPosition,subtreeTime);
                }
            }
        }
        else // only the one branch (topOfSlice) is in the slice now
        {
            lengths.push_back(topOfSlice.child->events.front()-topOfSlice.time);
            assert( lengths.back() > 0.0 );
        }
        length=accumulate(lengths.begin(),lengths.end(),0.0);
        double r = length * settings.fRand();

        numNodes=nodes.size();
        assert( numNodes == lengths.size() );
        assert( length > 0.0 );
        if (length < 1.0e-6)
        {
            cout << "Warning: total length of slice " << length << " reattaching at original position " << endl;
            return originalPosition;
        }
        if (debug)
            cout << "Num nodes: " << nodes.size() << " slice length " << length << endl;

        list<double>::iterator lengthIt = lengths.begin();
        for (list<BranchPoint>::iterator bpIt=nodes.begin(); bpIt != nodes.end(); ++bpIt)
        {
            r -= *lengthIt;
            if (r < 0.0) // then choose this reattachment position
            {
                r=-r; // how far into this segment to go
                BranchPoint bp=*bpIt;

                if (bp.isExistingBranchPoint)
                {
                    //bp.time = bp.child->node_event + *lengthIt - r;
                    bp.time = settings.nodeWidth - bp.child->node_event - *lengthIt + r;
                    if (bp.IsBefore(originalPosition) && bp.onPath)
                    {
                        topOfSlice = bp;
                    }
                    else
                    {
                        bp.time = bp.child->node_event + *lengthIt - r;
                        bp.child->node_event = bp.time;
                        bp.time = 0.0;
                    }
                }
                else
                {
                    double t0= bp.child==topOfSlice.child ? topOfSlice.time : bp.parent->time;
                    double t1= bp.child->events.empty() ? min(bp.child->time,subtreeTime) : bp.child->events.front();
                    assert( (t1-t0) == *lengthIt); // note this might fail because we're not taking subtreeTIme into account here
                    assert( r < *lengthIt );
                    bp.time=t0+r;
                    if ((bp.time > (originalPosition.isExistingBranchPoint ? originalPosition.child->time : originalPosition.time))  || !bp.onPath)
                    {
                        bp.child->events.clear();
                        bp.child->events.push_back(bp.time);
                    }
                    else
                    {
                        topOfSlice=bp;
                    }
                    if (bp.time < bp.parent->time || bp.time > bp.child->time)
                        throw 1;
                }

                return bp;
            }
            lengthIt++;
        }
        throw 1; // shouldn't get here!
    }



    // Choose a leaf at random to detach AND sample the Poisson events R,S,T
    TreeNode* GetRandomLeafAndSampleRST(BaseSettings &settings)
    {
        return root->GetRandomLeafAndSampleRST(settings, 0.0);
    }

    // Choose a subtree at random to detach AND sample the Poisson events R,S,T
    TreeNode* GetRandomSubtreeAndSampleRST(BaseSettings &settings)
    {
        return root->GetRandomSubtreeAndSampleRST(settings, 0.0);
    }

    // Sample Poisson events from the prior, used for inital tree building
    void SampleTPrior(BaseSettings &settings)
    {
        root->SampleT(settings, 0.0);
    }

    // Go through all the potential branch points calculating likelihoods
    // and then sample from the multinomial over them.
    void EvaluatePotentialBranchPoints(TreeNode* leaf,
                                       BaseSettings &s)
    {
        std::vector<BranchPoint> branchPoints; // potential branch points
        list<double> logProbs; // ...and their likelihoods
        root->EvaluatePotentialBranchPoints(branchPoints,logProbs,leaf,*zero, 0.0, 0.0, s);

//       cout << existing << " / " << probs.size() << endl;
        // sample from multinomial
        boost::numeric::ublas::vector<double> logProbsVector= listToVector(logProbs);
        boost::numeric::ublas::vector<double> probs = softmax(logProbsVector);
        int chosenBranchPointIndex = s.rdiscrete(probs);
        // attach at the choosen branch point
        branchPoints[chosenBranchPointIndex].AttachSubtree(leaf, s);
        // update numDescendants throughout the tree (could do this more efficiently!)
        zero->countLeaves(false);

        // the root may have changed
        root=*zero->children.begin();
    }

    /* Recurisvely propagate messages down the DDT.
    <param name="n">The node to propagate dowards from.</param>
    <param name="p">The parent of node n.</param> */
    void OutputTree(TreeNode &n, TreeNode &parent, ofstream &sw, int depth)
    {
        sw << depth << " "
           << n.Marg_Location.GetMean()[0] << " "
           <<  parent.Marg_Location.GetMean()[0] << " "
           <<   n.Marg_Location.GetMean()[1] << " "
           <<   parent.Marg_Location.GetMean()[1] << " "
           <<   n.time << " "
           <<   parent.time << endl;

        // Then recursively print the children.
        if (!n.isLeaf)
        {
            for (list<TreeNode*>::iterator i=n.children.begin(); i != n.children.end(); ++i)
                OutputTree(**i, n, sw, depth+1);
        }
    }

    /* Predict the density of the points in x_star */
    void OutputTree(const char* fname)
    {
        ofstream f;
        f.open(fname, ios::out);

        // Then send the messages recursively down the tree.
        OutputTree(*root, *zero, f, 0);
        f.close();
    }



    static double sqr(double x)
    {
        return x*x;
    }

    /*  Recurisvely propagate messages down the DDT.
     <param name="n">The node to propagate dowards from.</param>
     <param name="p">The parent of node n.</param> */
    boost::numeric::ublas::vector<double> Prediction(TreeNode &n, TreeNode &parent, list<boost::numeric::ublas::vector<double> > x_star, BaseSettings &settings)
    {
        // Here we calculate the probability of diverging on the edge leading to node n
        double s = parent.time;
        double t = n.time;
        double c = settings.c;
        double logProbNotDiverging, logProbDiverging;
        if (t == 0.0 || t <= s)
        {
            logProbDiverging = 0.0;
            logProbNotDiverging = log(0.0);
        }
        else
        {
            logProbNotDiverging = c * exp(lgamma(n.numDescendants - settings.alpha) - lgamma(n.numDescendants + 1 + settings.theta)) * (log((1.0 - t) / (1.0 - s)));
            logProbDiverging = log(1.0-exp(logProbNotDiverging));
            if (isnan(logProbDiverging))
                throw 1;
        }

        double divergenceTimes[settings.NumSamplesForImputation];
        double factor = c * exp(lgamma((double)n.numDescendants - settings.alpha) - lgamma((double)n.numDescendants + 1 + settings.theta));
        double Cs = 1.0 - pow(1.0 - s, factor);
        double Ct = 1.0 - pow(1.0 - t, factor);
        for (int i = 0; i < settings.NumSamplesForImputation; i++)
        {
            double u = Cs + settings.fRand() * (Ct - Cs); // draw uniform between Cs and Ct
            divergenceTimes[i] = 1.0 - pow(1.0 - u, 1.0 / factor);
            if (divergenceTimes[i] == 1.0)
                divergenceTimes[i] = (s + t) / 2.0;
        }
        double sigma = 1.0; // TODO: learning sigma?
        boost::numeric::ublas::vector<double> matrix[x_star.size()];
        for (int i=0; i<x_star.size(); i++)
        {
            matrix[i].resize(settings.NumSamplesForImputation);
        }
        for (int sample_index = 0; sample_index < settings.NumSamplesForImputation; sample_index++) // average over this
        {
            double divergenceTime = divergenceTimes[sample_index];
            boost::numeric::ublas::vector<double> divergenceMean = (parent.Marg_Location.GetMean() * (divergenceTime - s)
                    + n.Marg_Location.GetMean() * (t - divergenceTime)) / (t - s);
            boost::numeric::ublas::vector<double> divergenceVar = applyFunc((applyFunc(parent.Marg_Location.GetVariance(),sqrt) * (divergenceTime - s)
                    +applyFunc(n.Marg_Location.GetVariance(),sqrt) * (t - divergenceTime)) / (t - s),sqr);
            // Brownian bridge
            // divergenceVar += sigma * (divergenceTime - s) * (t - divergenceTime) / (t - s);
            divergenceVar += sigma * (1.0 - divergenceTime);
            int counter=0;
            for (list<boost::numeric::ublas::vector<double> >::iterator i=x_star.begin(); i != x_star.end(); ++i)
            {
                matrix[counter][sample_index] += sumVector(GaussianVector::GaussianLogProb(*i, divergenceMean, divergenceVar));
                assert(!isnan(matrix[counter][sample_index]));
            }
        }
        boost::numeric::ublas::vector<double> result(x_star.size());
        result &= 0.0;
        for (int i = 0; i < x_star.size(); i++)
        {
            result[i] += exp(logSumExp(matrix[i]) + logProbDiverging - log(settings.NumSamplesForImputation));
            assert(!isnan(result[i])); 
        }

        // Then recursively update the children.
        if (!n.isLeaf)
        {
            assert(exp(logProbNotDiverging) > 0.0 ); 
            boost::numeric::ublas::vector<double> divVar(settings.D);
            divVar &= sigma * (1.0 - n.time);
            double con = logProbNotDiverging - log(n.numDescendants + settings.theta);
            double probSplit = con + log(n.children.size() * settings.alpha + settings.theta);
            int counter=0;
            for (list<boost::numeric::ublas::vector<double> >::iterator i=x_star.begin(); i != x_star.end(); ++i)
            {
                result[counter] += exp(sumVector(GaussianVector::GaussianLogProb(*i, n.Marg_Location.GetMean(), divVar)) + probSplit);
                counter++;
            }
            for (list<TreeNode*>::iterator i=n.children.begin(); i != n.children.end(); ++i)
            {
                result += exp(con + log((*i)->numDescendants - settings.alpha)) * Prediction(**i, n, x_star, settings);
            }
        }


        return result;
    }


    /* Predict the density of the points in x_star. Vector returned has length equal to the number of these points,
     and is a _probability_ _not_ a _log prob_! */
    boost::numeric::ublas::vector<double> Prediction(list<boost::numeric::ublas::vector<double> > &x_star, BaseSettings &settings)
    {
        initialise(settings); // allocate memory for messages
        root->bpSweepUp(*zero); // belief propagation
        root->bpSweepDown(*zero);
        return Prediction(*root, *zero, x_star, settings);
    }

    boost::numeric::ublas::vector<double> PredictionSimple(list<boost::numeric::ublas::vector<double> > &x_star, BaseSettings &settings, int samples = 1000, int resampleEvery = 50)
    {

        boost::numeric::ublas::vector<double> logProbs[x_star.size()];
        boost::numeric::ublas::vector<double> result(x_star.size());
        result &= 0.0;
        for (int j=0; j<x_star.size(); j++)
            logProbs[j].resize(samples);

        for (int j=0; j<samples; j++)
        {
            if (j % resampleEvery == 0)
                InstantiateLocations(settings);
            GaussianVector g = root->PredictionSimple(settings);
            int counter=0;
            for (list<boost::numeric::ublas::vector<double> >::iterator i=x_star.begin(); i != x_star.end(); ++i)
            {
                logProbs[counter][j]=sumVector(g.GetLogProb(*i));
                counter++;
            }
        }

        for (int j=0; j<x_star.size(); j++)
            result[j]=logSumExp(logProbs[j])-log((double)samples);

        return result;
    }
};

#endif
