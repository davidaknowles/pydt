#define _USE_MATH_DEFINES

#ifndef TREENODE_H
#define TREENODE_H

#include <iostream>
#include <list>
#include <vector>
#include <map>
#include <cstdlib>
#include <math.h>
#include <boost/lexical_cast.hpp>
#include <boost/numeric/ublas/vector.hpp>
//#include <boost/thread.hpp>
#include <string>
#include <limits>
#include "GaussianVector.hpp"
#include "Settings.hpp"

#include <boost/math/special_functions/fpclassify.hpp>

using namespace std;

// Belief propagation message to "sample" from the Gaussian factor
// mean: incoming message from mean
// variance: known variance
GaussianVector SampleAverageConditional(GaussianVector mean, double variance)
{
    //cout << "Mean=" << mean << " prec=" << prec << endl;
    double prec = 1.0 / variance;
    if (mean.isPoint)
    {
        boost::numeric::ublas::vector<double> temp(mean.size());
        temp &= prec;
        return GaussianVector(mean.GetMean() * prec, temp);
    }
    else
    {
        boost::numeric::ublas::vector<double> R;
        R = (mean.Precision + prec) / prec;
        return GaussianVector(mean.MeanTimesPrecision / R, mean.Precision / R);
    }
};

class TreeNode;

// Potential branch point. If isExistingBranchPoint this this at an existing branch point (PYDT setting only)
// specified by "child". If not isExistingBranchPoint this is on the edge between "parent" and "child", at time
// "time"
class BranchPoint
{
public:
    TreeNode* parent;
    TreeNode* child;
    double time;
    bool isExistingBranchPoint;
    bool onPath;

    void AttachSubtree(TreeNode* subtree, BaseSettings &settings);

    string ToString();

    bool IsBefore(BranchPoint &other);

    BranchPoint() {};

    BranchPoint(TreeNode* parent, TreeNode* child, double time, bool isExistingBranchPoint, bool onPath) :
        parent(parent),
        child(child),
        time(time),
        isExistingBranchPoint(isExistingBranchPoint),
        onPath(onPath) {};

    BranchPoint(TreeNode* parent, TreeNode* child, double time, bool isExistingBranchPoint) :
        parent(parent),
        child(child),
        time(time),
        isExistingBranchPoint(isExistingBranchPoint) {};


};

template <typename T>
T listIndexOf(list<T> l, int index)
{
    int counter=0;
    for (typename list<T>::iterator i=l.begin(); i != l.end(); ++i)
    {
        if (counter==index)
            return *i;
        counter++;
    }
    throw 1;
}

// Return whether the list has N elements
template <typename T>
bool ListHasOnlyNElements(list<T> l, int N)
{
    //typename T;
    typename list<T>::iterator it;
    it=l.begin();
    for (int i=0; i<N; i++)
    {
        if (it==l.end())
            return false; // list has i elements
        it++;
    }
    if (it!=l.end())
        return false; // list has more than 1 element
    return true;
}

// Return whether the list has only one element
template <typename T>
bool ListHasOnlyOneElement(list<T> l)
{
    return ListHasOnlyNElements<T>(l, 1);
}

// Class representing a node in the tree
class TreeNode {

    friend ostream &operator<<(ostream &output, const TreeNode &node)
    {
        output << node.time << endl;
        return output;
    }

public:

    // Divergence times
    double time;

    // List of children
    list<TreeNode*> children;

  // Ideally we wouldn't store this but makes some operations much easier
  TreeNode* parent; 

    // Number of leaf nodes found down the tree from here
    int numDescendants;

    // Whether this is a leaf node
    bool isLeaf;

    // A text label for this node
    string label;

    // Marginal location, also used for instantiated location
    GaussianVector Marg_Location;

    // BP message from this to parent
    GaussianVector Msg_Normal_ParentLocation;

    // BP message from parent to this
    GaussianVector  Msg_Normal_Location;

    // list of Poisson events on the branch from parent to this
    list<double> events;

    // list of Poisson events on the branch from parent to this
    double node_event;

    // Constructor: t is time, il = isLeaf
    TreeNode (double t, bool il) {
        time = t;
        isLeaf=il;
        numDescendants = 1;
        label="";
        node_event = 0.0;
    };

  void SetupParents(TreeNode* theparent, bool check=false){
    if (check)
      assert(parent==theparent); 
    parent=theparent; 
    for (list<TreeNode*>::iterator i=children.begin(); i != children.end(); ++i)
      (*i)->SetupParents(this,check); 
  }

    // Check that the BP messages are consistent with the marginal
    void CheckConsistency()
    {
        GaussianVector temp = Msg_Normal_Location;
        if (!isLeaf)
        {
            int childCount = 0;
            for (list<TreeNode*>::iterator i=children.begin(); i != children.end(); ++i)
            {
                temp *= (*i)->Msg_Normal_ParentLocation;
                (*i)->CheckConsistency();
                childCount++;
            }
            if (childCount <= 1)
                throw 1;
            cout << temp << "=" << Marg_Location << endl;
        }

    }

    // Delete a subtree recursively free-ing associated memory
    static void deleteSubtree(TreeNode* subtree)
    {
        for (list<TreeNode*>::iterator i=subtree->children.begin(); i != subtree->children.end(); ++i)
        {
            deleteSubtree(*i);
        }
        delete subtree;
    };

    // Get a deep copy of this subtree recursively
    TreeNode* DeepCopy()
    {
        TreeNode* result = new TreeNode(time, isLeaf);
        result->numDescendants=numDescendants;
        result->label=label;
        if (isLeaf)
            result->Marg_Location=Marg_Location; // observed data!
        // NB: not copying messages!
        for (list<TreeNode*>::iterator i=children.begin(); i != children.end(); ++i)
        {
            result->children.push_back((*i)->DeepCopy());
        }
        return result;
    }

    // Used for aggregating depth statistics about the tree: how many leaves are there
    // at each depth
    void LeafDepthHist(map<int,int> &hist, int depth = 0)
    {
        if (!isLeaf)
        {
            for (list<TreeNode*>::iterator i=children.begin(); i != children.end(); ++i)
            {
                (*i)->LeafDepthHist(hist,depth+1);
            }
        }
        else
        {
            if (hist.count(depth)==0)
                hist[depth]=0;
            hist[depth]++;
        }
    }

    int sumBranchingFactor()
    {
        if (!isLeaf)
        {
            int sum=children.size(); 
            for (list<TreeNode*>::iterator i=children.begin(); i != children.end(); ++i)
            {
                sum += (*i)->sumBranchingFactor();
            }
            return sum; 
        }
        else
        {
            return 0; 
        }
    }

    int maxBranchingFactor()
    {
        if (!isLeaf)
        {
            int maxbf=children.size(); 
            for (list<TreeNode*>::iterator i=children.begin(); i != children.end(); ++i)
            {
                maxbf = max( (*i)->maxBranchingFactor(), maxbf ); 
            }
            return maxbf; 
        }
        else
        {
            return 0; 
        }
    }

    // Detach a subtree. Returns who the root should be. logProb is the (hypothetical) probability of attaching here
    TreeNode* DetachSubtree(TreeNode* subtree, TreeNode* parent, BranchPoint* originalPosition = NULL)
    {
        for (list<TreeNode*>::iterator i=children.begin(); i != children.end(); ++i)
        {
            if (subtree==*i) // this child is the subtree to be removed
            {
                if (originalPosition != NULL) // record where the subtree was
                {
                    originalPosition->parent=parent;
                    originalPosition->child=this;
                    originalPosition->time=(*i)->time;
                    originalPosition->isExistingBranchPoint=true;
                }
                children.remove(*i); // remove it from our list of children
                numDescendants -= subtree->numDescendants;
                // if removing this child meant we now only have one child
                // then delete this
                if (ListHasOnlyOneElement<TreeNode*>(children))
                {
                    if (originalPosition != NULL)
                    {
                        originalPosition->parent=parent;
                        originalPosition->child=children.front();
                        originalPosition->time=time;
                        originalPosition->isExistingBranchPoint=false;
                    }
                    TreeNode* toreturn= children.front();
		    toreturn->parent=parent; 
                    delete this;
                    return toreturn;
                }
                else
                {
                    return this;
                }
            }
            TreeNode* temp = (*i)->DetachSubtree(subtree,this,originalPosition);
            if (temp!=NULL) // if subtree was detached below this child...
            {
                numDescendants -= subtree->numDescendants;
                (*i) = temp; // the child may have been deleted
                return this;
            }
        }
        return NULL; // the subtree to detach was not in the subtree rooted at this node
    }

    // Initialise messages before performing a sweep of BP
    void initialise(BaseSettings &settings)
    {
        Msg_Normal_ParentLocation.SetToUniform(settings.D);
        Msg_Normal_Location.SetToUniform(settings.D);
        if (!isLeaf)
        {
            Marg_Location.SetToUniform(settings.D);
            for (list<TreeNode*>::iterator i=children.begin(); i != children.end(); ++i)
                (*i)->initialise(settings);
        }
    };

    // Add a leaf according to the generative process. Returns the attachment depth
    int AddChild(BaseSettings& s, TreeNode* leaf, int depth=0) {
        double u = s.fRand(0,numDescendants+s.theta);
        u -= s.theta+s.alpha*(double)(int)children.size();
        numDescendants++;
        if (u < 0) // create a new branch at an existing branch point (PYDT only)
        {
            children.push_back(leaf);
            return depth;
        }
        else
        {
            for (list<TreeNode*>::iterator i=children.begin(); i != children.end(); ++i)
            {
                u -= (*i)->numDescendants - s.alpha;
                if (u < 0) // go down this branch
                {
                    // decide whether to diverge on this branch, and if so, when
                    double At = s.A(time)+exp(s.logDivergenceRateFactor((*i)->numDescendants))*(-log(s.fRand(0,1)));
                    double td= s.invA(At);
                    if (td>(*i)->time) // don't diverge off this branch
                        return (*i)->AddChild(s,leaf,depth+1);
                    else // diverge
                    {
                        // create new internal node
                        TreeNode * temp = new TreeNode(td,false);
                        temp->children.push_back(*i);
                        temp->children.push_back(leaf);
                        temp->numDescendants = (*i)->numDescendants+1;
                        *i=temp;
                        return depth+1;
                    }
                }
            }
        }
        throw 1; // should never get here!
    };

    // NOTE: this assumes sigma=1
    GaussianVector PredictionSimple(BaseSettings& s) {
        double u = s.fRand(0,numDescendants+s.theta);
        u -= s.theta+s.alpha*(double)(int)children.size();
        GaussianVector result;
        if (u < 0) // create a new branch at an existing branch point (PYDT only)
        {
            //result= new GaussianVector();
            boost::numeric::ublas::vector<double> variance(s.D);
            variance &= 1.0 - time;
            result.SetMeanAndVariance(Marg_Location.GetMean(), variance);
            return result;
        }
        else
        {
            for (list<TreeNode*>::iterator i=children.begin(); i != children.end(); ++i)
            {
                u -= (*i)->numDescendants - s.alpha;
                if (u < 0) // go down this branch
                {
                    // decide whether to diverge on this branch, and if so, when
                    double At = s.A(time)+exp(s.logDivergenceRateFactor((*i)->numDescendants))*(-log(s.fRand(0,1)));
                    double td= s.invA(At);
                    if (td>(*i)->time) // don't diverge off this branch
                        return (*i)->PredictionSimple(s);
                    else // diverge
                    {
                        // create new internal node
                        //result= new GaussianVector();
                        boost::numeric::ublas::vector<double> variance(s.D);
                        variance &= 1.0 - td;
                        boost::numeric::ublas::vector<double> mean = (Marg_Location.GetMean()*(td-time)) + ((*i)->Marg_Location.GetMean()*((*i)->time-td));
                        mean /= (*i)->time - time;
                        result.SetMeanAndVariance(mean, variance);
                        return result;
                    }
                }
            }
        }
        throw 1; // should never get here!
    };

    // Try to attach the subtree using the generative process
    // Return whether this was successful
    // isZero: whether this is the very first node in the tree (above the root) so branching
    // directly from here is not allowed
    // Returns whether we successfully attached the subtree
  bool AttachSubtree(BaseSettings& s, TreeNode* subtree,BranchPoint &where, bool isZero = true) {

        if (subtree->time < time)
            return false;
        double u = s.fRand(0.0,numDescendants+s.theta) ;
        u -= s.theta+s.alpha*(int)children.size();
        // form a new branch from this existing branch point
        if ((!isZero) && u < 0.0)
        {
            children.push_back(subtree);
	    where.parent=this;
	    where.child=NULL;
	    where.time=time; 
	    where.isExistingBranchPoint=true;
            numDescendants+=subtree->numDescendants;
            //cout << "new branch" << endl;
            return true;
        }
        else
        {
            for (list<TreeNode*>::iterator i=children.begin(); i != children.end(); ++i)
            {
                u -= (*i)->numDescendants-s.alpha;
                if (isZero || u < 0.0) // go down this branch
                {
                    // effective modification to the rate as a result of the multiple datapoints in the subtree
                    double logDivRateFactor=s.logDivergenceRateFactor((*i)->numDescendants);
                    double At = s.A(time)+exp(logDivRateFactor)*(-log(s.fRand()));
                    double td= s.invA(At);
                    if (td>(*i)->time) 	    {

		      bool temp=(*i)->AttachSubtree(s,subtree,where,false);
		      if (temp) // if the subtree was attached below us in the tree then need to update numDescendants
                        {
			  numDescendants+=subtree->numDescendants;
                        }
                        return temp;
                    }
                    else
                    {
                        if (td<subtree->time) // check that the branch length will be positive
                        {
                            //cout << "attached at " << td << endl;
                            TreeNode * temp = new TreeNode(td,false);
                            temp->children.push_back(*i);
			    where.child=*i; 
                            temp->children.push_back(subtree);
                            temp->numDescendants = (*i)->numDescendants+subtree->numDescendants;
                            children.remove(*i);
                            children.push_back(temp);
                            numDescendants+=subtree->numDescendants;
			    where.parent=this;
			    where.time=td;
			    where.isExistingBranchPoint=false;
                            return true;
                        }
                        else
                            return false; // did not attach
                    }
                }
            }
        }
        throw 1;
    };

    // Get the probability of having attached subtree at its current location
    bool GetProbOfAttachSubtree(BaseSettings& s, TreeNode* subtree, double &logProb)
    {
        // If subtree is one of our children find the factor for this
        for (list<TreeNode*>::iterator i=children.begin(); i != children.end(); ++i)
        {
            if ((*i)==subtree)
            {
                if (ListHasOnlyNElements(children,2)) // adding subtree resulted in the creation of this node
                    logProb += log(s.a(time)) + s.singleTerm(numDescendants-subtree->numDescendants); // prob of diverging here
                else
                    logProb += log(s.theta+s.alpha*(int)(children.size()-1)) - log(numDescendants-subtree->numDescendants+s.theta);
                return true;
            }
        }


        for (list<TreeNode*>::iterator i=children.begin(); i != children.end(); ++i)
        {
            double temp = 0.0;
            bool attached = (*i)->GetProbOfAttachSubtree(s, subtree, temp);
            if (attached)
            {
                //cout << logProb << endl;
                logProb += (s.A(time)-s.A((*i)->time)) * exp(s.singleTerm((*i)->numDescendants-subtree->numDescendants)); // prob not diverging until td
                // cout << logProb << endl;
                logProb += log((*i)->numDescendants-subtree->numDescendants-s.alpha) - log(numDescendants-subtree->numDescendants+s.theta); // prob going down this branch
                //cout << logProb << endl;
                logProb += temp;
                // cout << logProb << endl;
                return true;
            }
        }

        return false;
    };

    // Output in Newick format, with divergence times
    string newick(double parentTime)
    {
        string res; // result string
        if (isLeaf) // just put the label if this is a leaf
        {
            res = label;
        }
        else // build up the string
        {
            res = "(";
            for (list<TreeNode*>::iterator i=children.begin(); i != children.end(); )
            {
                res += (*i)->newick(time);
                ++i;
                if (i != children.end())
                    res += ",";
            }
            res += ")";
        }
        // append the divergence time
        // note that some plotting packages expect time-parentTime here.
        res += ":" + boost::lexical_cast<string>(time) + "-" + boost::lexical_cast<string>(parentTime);

        return res;
    };

    string newick_struct(BaseSettings &settings)
    {
        string res; // result string
        if (isLeaf) // just put the label if this is a leaf
        {
            res = label;
        }
        else // build up the string
        {
            res = "(";
            for (list<TreeNode*>::iterator i=children.begin(); i != children.end(); )
            {
                res += (*i)->newick_struct(settings);
                ++i;
                if (i != children.end())
                    res += ",";
            }
            res += ")";
        }
        res += ":" + boost::lexical_cast<string>(LocalEvidenceStructure(settings)); 

        return res;
    };

  string newick_times_evidence(double parenttime, BaseSettings &settings)
    {
        string res; // result string
        if (isLeaf) // just put the label if this is a leaf
        {
            res = label;
        }
        else // build up the string
        {
            res = "(";
            for (list<TreeNode*>::iterator i=children.begin(); i != children.end(); )
            {
	      res += (*i)->newick_times_evidence(time,settings);
                ++i;
                if (i != children.end())
                    res += ",";
            }
            res += ")";
        }
        res += ":" + boost::lexical_cast<string>(LocalEvidenceTimes(parenttime,settings)); 

        return res;
    };


    // Output tree in newick format (structure only)
    string newick()
    {
        string res;
        if (isLeaf)
        {
            res = label;
        }
        else
        {
            res = "(";
            for (list<TreeNode*>::iterator i=children.begin(); i != children.end(); )
            {
                res += (*i)->newick();
                ++i;
                if (i != children.end())
                    res += ",";
            }
            res += ")";
        }
        //if (!events.empty())
        //res += ":" + boost::lexical_cast<string>(events.front());

        return res;
    };

    //  Output tree in newick format including divergence times AND locations
    string newick2(double parentTime)
    {
        string res;
        if (isLeaf)
        {
            res = label;
        }
        else
        {
            res = "(";
            for (list<TreeNode*>::iterator i=children.begin(); i != children.end(); )
            {
                res += (*i)->newick2(time);
                ++i;
                if (i != children.end())
                    res += ",";
            }
            res += ")";
        }
        stringstream s;
        s << Marg_Location.GetMean();
        res += ":" + s.str() + "-" + boost::lexical_cast<string>(time);
        return res;
    };

    // Count the number of leaves.
    // check: whether to throw an exception if numDescendants is incorrect
    int countLeaves(bool check = false, bool countAll = false)
    {
        if (isLeaf)
            return 1;
        int leaves = 0;
        int childrenCount = 0;
        for (list<TreeNode*>::iterator i=children.begin(); i != children.end(); ++i,childrenCount++)
            leaves += (*i)->countLeaves(check, countAll);
        if (check && (leaves != numDescendants))
        {
            cerr << "Number of leaves: " << leaves << " cached number=" << numDescendants << endl;
	    cout << newick() << endl; 
            throw 1;
        }
        else if (check && childrenCount==1)
            throw 1;
        else if (!countAll)
            numDescendants=leaves;

        return leaves + (countAll ? 1 : 0);
    };




    // Upwards sweep of sampling
    void sampleSweepUp(TreeNode& parent, BaseSettings &settings)
    {
        parent.Marg_Location /= Msg_Normal_ParentLocation; // remove the current contribution
        if (!Marg_Location.isPoint)
            Marg_Location.SetPoint(Marg_Location.Sample(settings.gen)); // sample point
        if (isLeaf)
        {
            // calculate the new message
            Msg_Normal_ParentLocation = SampleAverageConditional(Marg_Location, time-parent.time);
        }
        else
        {
            for (list<TreeNode*>::iterator i=children.begin(); i != children.end(); ++i)
                (*i)->sampleSweepUp(*this, settings); // recurse down the structure
            // ... then calculate messages
            Msg_Normal_ParentLocation = SampleAverageConditional(Marg_Location, time-parent.time);
        }
        // update the parent marginal
        parent.Marg_Location *= Msg_Normal_ParentLocation;
        // cout << Msg_Normal_ParentLocation << time - parent.time << endl;
    };


    // Upwards sweep of belief propagation
    void bpSweepUp(TreeNode& parent)
    {
        parent.Marg_Location /= Msg_Normal_ParentLocation; // remove the current contribution
        if (isLeaf)
        {
            // calculate the new message
            Msg_Normal_ParentLocation = SampleAverageConditional(Marg_Location, time-parent.time);
        }
        else
        {
            for (list<TreeNode*>::iterator i=children.begin(); i != children.end(); ++i)
                (*i)->bpSweepUp(*this); // recurse down the structure
            // ... then calculate messages
            Msg_Normal_ParentLocation = SampleAverageConditional(Marg_Location / Msg_Normal_Location, time-parent.time);
        }
        // update the parent marginal
        parent.Marg_Location *= Msg_Normal_ParentLocation;
        // cout << Msg_Normal_ParentLocation << time - parent.time << endl;
        if (isnan(Msg_Normal_ParentLocation.Precision[0])) throw 1;
    };

    // Upwards sweep of belief propagation starting at subtree attachment or detachment
    void bpSweepUp(int numDescendantsToAdd)
    {
      if (parent!=NULL){
        parent->Marg_Location /= Msg_Normal_ParentLocation; // remove the current contribution
	// calculate the new message
	Msg_Normal_ParentLocation = SampleAverageConditional(isLeaf ? Marg_Location : (Marg_Location / Msg_Normal_Location), time-parent->time);
	// update the parent marginal
        parent->Marg_Location *= Msg_Normal_ParentLocation;
	parent->numDescendants += numDescendantsToAdd; 
        assert(!isnan(Msg_Normal_ParentLocation.Precision[0])); 
	parent->bpSweepUp(numDescendantsToAdd); 
      }
    };


    // Upwards sweep of belief propagation
  double bpSweepUp2(GaussianVector &msg_up, bool isRoot,int D)
  {
    if (isLeaf) {
      msg_up = Marg_Location;
      return 0.0; 
    }
    else
      {
	double ml=0.0;
	boost::numeric::ublas::vector<double> prod_v(D), sum_m2_over_v(D);
	prod_v &= 1.0 ;
	sum_m2_over_v &= 0.0; 
	msg_up.SetToUniform(D); 
	for (list<TreeNode*>::iterator i=children.begin(); i != children.end(); ++i){
	  GaussianVector up_msg; 
	  ml += (*i)->bpSweepUp2(up_msg,false,D); // recurse down the structure
	  GaussianVector norm_up=SampleAverageConditional(up_msg, (*i)->time - time); 
	  msg_up *= norm_up; 
	  boost::numeric::ublas::vector<double> mi=norm_up.GetMean(); 
	  boost::numeric::ublas::vector<double> ti=norm_up.GetVariance(); 
	  prod_v *= ti; 
	  sum_m2_over_v += mi*mi/ti; 
	}
	if (isRoot){
	  GaussianVector prior;
	  prior.SetToUniform(D); 
	  prior.Precision &= 1.0 / time; 
	  msg_up *= prior; 
	  boost::numeric::ublas::vector<double> mi=prior.GetMean(); 
	  boost::numeric::ublas::vector<double> ti=prior.GetVariance(); 
	  prod_v *= ti; 
	  sum_m2_over_v += mi*mi/ti; 
	}
	boost::numeric::ublas::vector<double> m =msg_up.GetMean(); 
	boost::numeric::ublas::vector<double> v=msg_up.GetVariance(); 
	ml += -((double)children.size()- (isRoot?0.0:1.0))*(double)D*GaussianVector::lnSqrt2Pi + .5 * ( sumVector(applyFunc(v,log)) - sumVector(applyFunc(prod_v,log)) + sumVector(m*m/v) - sumVector(sum_m2_over_v) ); 
	if (isnan(ml)) throw 1; 
	return ml; 
      }
  };
  
  double local_factor(bool excludeLast,int D,TreeNode* child_to_exclude=NULL){
    if (isLeaf) throw 1;
    boost::numeric::ublas::vector<double> prod_v(D), sum_m2_over_v(D);
    prod_v &= 1.0 ;
    sum_m2_over_v &= 0.0; 
    GaussianVector msg_up; 
    msg_up.SetToUniform(D); 
    int counter=0,nchildren=children.size();
    bool found_child_to_exclude=false; 
    for (list<TreeNode*>::iterator i=children.begin(); i != children.end(); ++i){
      if (excludeLast && counter==nchildren-1)
	continue; 
      if (child_to_exclude!=NULL && child_to_exclude==*i){
	found_child_to_exclude=true;
	continue;
      }
      GaussianVector norm_up=(*i)->Msg_Normal_ParentLocation; 
      msg_up *= norm_up; 
      boost::numeric::ublas::vector<double> mi=norm_up.GetMean(); 
      boost::numeric::ublas::vector<double> ti=norm_up.GetVariance(); 
      prod_v *= ti; 
      sum_m2_over_v += mi*mi/ti; 
      counter++;
    }
    if (child_to_exclude)
      assert(found_child_to_exclude); 
    msg_up *= Msg_Normal_Location;
    boost::numeric::ublas::vector<double> mi=Msg_Normal_Location.GetMean(); 
    boost::numeric::ublas::vector<double> ti=Msg_Normal_Location.GetVariance(); 
    prod_v *= ti; 
    sum_m2_over_v += mi*mi/ti; 
    boost::numeric::ublas::vector<double> m =msg_up.GetMean(); 
    boost::numeric::ublas::vector<double> v=msg_up.GetVariance(); 
    double ml= -(double)counter*(double)D*GaussianVector::lnSqrt2Pi + .5 * ( sumVector(applyFunc(v,log)) - sumVector(applyFunc(prod_v,log)) + sumVector(m*m/v) - sumVector(sum_m2_over_v) ); 
    assert(!isnan(ml)); 
    return ml;

  }

  static double local_factor(list<GaussianVector*> &in,int D){
    boost::numeric::ublas::vector<double> prod_v(D), sum_m2_over_v(D);
    prod_v &= 1.0 ;
    sum_m2_over_v &= 0.0; 
    GaussianVector msg_up; 
    msg_up.SetToUniform(D); 
    for (list<GaussianVector*>::iterator i=in.begin(); i != in.end(); ++i){
      GaussianVector* norm_up = *i; 
      msg_up *= *norm_up; 
      boost::numeric::ublas::vector<double> mi=norm_up->GetMean(); 
      boost::numeric::ublas::vector<double> ti=norm_up->GetVariance(); 
      prod_v *= ti; 
      sum_m2_over_v += mi*mi/ti; 
    }
    boost::numeric::ublas::vector<double> m =msg_up.GetMean(); 
    boost::numeric::ublas::vector<double> v=msg_up.GetVariance(); 
    return -((double)in.size()-1.0)*(double)D*GaussianVector::lnSqrt2Pi + .5 * ( sumVector(applyFunc(v,log)) - sumVector(applyFunc(prod_v,log)) + sumVector(m*m/v) - sumVector(sum_m2_over_v) ); 

  }



    // Downward sweep of belief propagation
    void bpSweepDown(TreeNode& parent)
    {
        if (isLeaf)
        {
            Msg_Normal_Location = SampleAverageConditional(parent.Marg_Location / Msg_Normal_ParentLocation, time-parent.time);
        }
        else
        {
            Marg_Location /= Msg_Normal_Location; // remove the current contribution
            // calculate the ne message
            Msg_Normal_Location = SampleAverageConditional(parent.Marg_Location / Msg_Normal_ParentLocation, time-parent.time);
            Marg_Location *= Msg_Normal_Location; // add in new contribution
            for (list<TreeNode*>::iterator i=children.begin(); i != children.end(); ++i)
                (*i)->bpSweepDown(*this); // recurse down the structure
            //cout << "Marginal: " << Marg_Location << endl;
        }
        if (isnan(Msg_Normal_Location.Precision[0])) throw 1;
        //cout << "Marginal: " << Marg_Location << endl;
    };

    // Sample synthetic data at the leaves
    void sampleData(boost::numeric::ublas::vector<double> &parentLocation, double parentTime, BaseSettings &s)
    {
        // sample the diffusion process for this branch
        boost::numeric::ublas::vector<double> ones(parentLocation.size());
        ones &= 1.0;
        boost::numeric::ublas::vector<double> temp = GaussianVector::Sample(s.gen, parentLocation, ones * (time-parentTime));
        if (isLeaf)
        {
            Marg_Location.SetPoint(temp);
            //cout << label << " " << temp << endl;
        }
        else
        {
            for (list<TreeNode*>::iterator i=children.begin(); i != children.end(); ++i)
                (*i)->sampleData(temp, time, s);
        }
    };

  double LocalEvidenceTimes(double parent_time, BaseSettings &settings){
    double lnEvidence = 0.0;
    if (!isLeaf)
      {
	lnEvidence += (settings.A(parent_time) - settings.A(time))*settings.H(numDescendants-1);
	lnEvidence += log(settings.a(time));
      }
    assert(!isnan(lnEvidence));
    return lnEvidence;
  }

  double LocalEvidenceStructure(BaseSettings &settings){
    double lnEvidence = 0.0;
    if (!isLeaf)
      {
	// is this correct?
	int l=1;
	for (list<TreeNode*>::iterator i=children.begin(); i != children.end(); ++i)
	  {
	    if (l >= 3)
	      lnEvidence += log(settings.theta+((double)l-1.0)*settings.alpha);
	    lnEvidence += lgamma((*i)->numDescendants-settings.alpha);
	    l++;
	  }
	lnEvidence -= lgamma(numDescendants+settings.theta);
	int numBranches = l - 1; 
	lnEvidence -= ((double)numBranches-1.0)*lgamma(1.0-settings.alpha);
      }
    assert(!isnan(lnEvidence));
    return lnEvidence;
  }

    // log evidence contribution for the structure and divergence times
    double LogEvidenceStructure(double parent_time, BaseSettings &settings)
    {
      double lnEvidence = LocalEvidenceStructure(settings); 
      lnEvidence += LocalEvidenceTimes(parent_time,settings); 
      for (list<TreeNode*>::iterator i=children.begin(); i != children.end(); ++i)
	lnEvidence += (*i)->LogEvidenceStructure(time, settings);
      return lnEvidence;
    }

  // additional_descendants is used to simulated the effect of adding a subtree
  double LogEvidenceStructureUp(BaseSettings &settings,int additional_descendants)
    {
      double lnEvidence = 0.0; 
      if (parent != NULL) { // i.e. stop at zero!
	numDescendants += additional_descendants; 
	lnEvidence += LocalEvidenceStructure(settings); 
	lnEvidence += LocalEvidenceTimes(parent->time,settings); 
	lnEvidence += parent->LogEvidenceStructureUp(settings,additional_descendants);
	numDescendants -= additional_descendants; 
      }
      return lnEvidence;
    }

  double LogEvidenceStructureOnlyUp(BaseSettings &settings,int additional_descendants)
    {
      double lnEvidence = 0.0; 
      if (parent != NULL) { // i.e. stop at zero!
	numDescendants += additional_descendants; 
	cout << "Up: " << LocalEvidenceStructure(settings) << "add: " << additional_descendants << endl; 
	lnEvidence += LocalEvidenceStructure(settings); 
	lnEvidence += parent->LogEvidenceStructureOnlyUp(settings,additional_descendants);
	numDescendants -= additional_descendants; 
      }
      return lnEvidence;
    }

  double LogEvidenceTimesOnlyUp(BaseSettings &settings,int additional_descendants)
    {
      double lnEvidence = 0.0; 
      if (parent != NULL) { // i.e. stop at zero!
	numDescendants += additional_descendants; 
	lnEvidence += LocalEvidenceTimes(parent->time,settings); 
	cout << "Time up: " << LocalEvidenceTimes(parent->time,settings) << " add " << additional_descendants << endl; 
	lnEvidence += parent->LogEvidenceTimesOnlyUp(settings,additional_descendants);
	numDescendants -= additional_descendants; 
      }
      return lnEvidence;
    }


    double LogEvidenceTimes(double parent_time, BaseSettings &settings)
    {
      double lnEvidence = LocalEvidenceTimes(parent_time,settings); 
      cout << "time: " << time << " parent: " << parent_time << " ev: " << lnEvidence << endl; 
      if (!isLeaf)
	for (list<TreeNode*>::iterator i=children.begin(); i != children.end(); ++i)
	  lnEvidence += (*i)->LogEvidenceTimes(time, settings);
      return lnEvidence;
    }

    double LogEvidenceStructureOnly(BaseSettings &settings)
    {
      double lnEvidence = LocalEvidenceStructure(settings);
      if (!isLeaf)
	for (list<TreeNode*>::iterator i=children.begin(); i != children.end(); ++i)
	  lnEvidence += (*i)->LogEvidenceStructureOnly(settings);
      return lnEvidence;
    }

    // log evidence contribution of the locations, calculated EP style (probably can do this more straightforwardly!)
    // I've compared this to Infer.NET to make sure it's correct
    double LogEvidence(TreeNode& parent)
    {
        if (isLeaf)
        {
            if (!Marg_Location.isPoint)
                throw 1; // not implemented yet!
            double result = sumVector(Msg_Normal_Location.GetLogProb(Marg_Location.GetMean()));
            if (isnan(result))
                throw 1;
            return result;
        }

        double lnEvidence = (1.0-(double)children.size())*Marg_Location.GetLogNormalizer() - Msg_Normal_Location.GetLogNormalizer();
        //std::cout << "Marg_Location " << Marg_Location << std::endl;
        for (list<TreeNode*>::iterator i=children.begin(); i != children.end(); ++i)
        {
            GaussianVector mean = Marg_Location / (*i)->Msg_Normal_ParentLocation;
            lnEvidence += mean.GetLogNormalizer();
            //std::cout << "mean " << Marg_Location << std::endl;
            //std::cout << "Msg_Normal_ParentLocation " << (*i)->Msg_Normal_ParentLocation << std::endl;
            if (isnan(lnEvidence))
                throw 1;
            lnEvidence += (*i)->LogEvidence(*this);
        }
        if (isnan(lnEvidence))
            throw 1;
        return lnEvidence;
    }


    // log evidence contribution of the locations, calculated EP style (probably can do this more straightforwardly!)
    // I've compared this to Infer.NET to make sure it's correct
    double LogEvidence2(TreeNode& parent)
    {
        if (isLeaf)
        {
            if (!Marg_Location.isPoint)
                throw 1; // not implemented yet!
            return 0.0; 
        }

        double lnEvidence = (1.0-(double)children.size())*Marg_Location.GetLogNormalizer() - Msg_Normal_Location.GetLogNormalizer();
        //std::cout << "Marg_Location " << Marg_Location << std::endl;
        for (list<TreeNode*>::iterator i=children.begin(); i != children.end(); ++i)
        {
            GaussianVector mean = Marg_Location / (*i)->Msg_Normal_ParentLocation;
            lnEvidence += mean.GetLogNormalizer();
            //std::cout << "mean " << Marg_Location << std::endl;
            //std::cout << "Msg_Normal_ParentLocation " << (*i)->Msg_Normal_ParentLocation << std::endl;
            if (isnan(lnEvidence))
                throw 1;
            lnEvidence += (*i)->LogEvidence(*this);
        }
        if (isnan(lnEvidence))
            throw 1;
        return lnEvidence;
    }


    // Get a random subtree of the tree rooted at this node
    TreeNode* GetRandomSubtree(BaseSettings &settings, bool isRoot=true)
    {
        // if we got to a leaf then we must choose this
        if (isLeaf)
            return this;
        // return this subtree?
        double probMe = settings.subtreeInflationFactor / (settings.subtreeInflationFactor + (double)numDescendants);
        if ((settings.fRand() < probMe) && !isRoot) // the subtree rooted at this node will be detached
            return this;

        // if not, randomly select a branch to go down
        int temp=settings.iRand(0, numDescendants-1);
        // choose a branch at random, with probability proportional to the number of leaves
        // under that branch
        int counter = 0;
        for (list<TreeNode*>::iterator i=children.begin(); i != children.end(); ++i)
        {
            counter+=(*i)->numDescendants;
            if (counter>=temp)
                return (*i)->GetRandomSubtree(settings,false);
        }
        throw 1; // shouldn't ever get here
    }

    void GetAllNodes(list<TreeNode*> &allNodes, bool isRoot = true)
    {
        if (!isRoot)
            allNodes.push_back(this);
        for (list<TreeNode*>::iterator i=children.begin(); i != children.end(); ++i)
        {
            (*i)->GetAllNodes(allNodes,false);
        }
    }

    bool GetBranchesInSlice(BaseSettings &settings, list<BranchPoint> &branches, list<double> &lengths, TreeNode& parent, BranchPoint& originalPosition, double subtreeTime)
    {
        if (!events.empty()) // there is a previously rejected reattachment position on this branch, so the slice ends here
        {
            branches.push_back(BranchPoint(&parent,this,0.0,false,originalPosition.child==this));
            lengths.push_back((*events.begin())-parent.time);
            assert( lengths.back() > 0.0 );
            return originalPosition.child==this;
        }
        else if (node_event == 0.0) // there has NOT been a rejected reattachment position at this node 
        {
            bool onPath=originalPosition.child==this;

            for (list<TreeNode*>::iterator i=children.begin(); i != children.end(); ++i)
            {
                if (time < subtreeTime)
                {
                    if ((*i)->GetBranchesInSlice(settings, branches,lengths,*this,originalPosition,subtreeTime))
                        onPath=true;
                }
                if (*i == originalPosition.child)
                    onPath=true;
            }
            if (subtreeTime > parent.time)
            {
                branches.push_back(BranchPoint(&parent,this,0.0,false,onPath));
                lengths.push_back(min(time,subtreeTime)-parent.time);
                assert( lengths.back() > 0.0 );
            }

            if (settings.multifurcating() && subtreeTime > time)
            {
                branches.push_back(BranchPoint(&parent,this,0.0,true,onPath));
                lengths.push_back(settings.nodeWidth);
            }

            return onPath;
        }
        else
        {
            bool onPath=originalPosition.child==this;

            if (settings.multifurcating() && subtreeTime > time)
            {
                branches.push_back(BranchPoint(&parent,this,0.0,true,onPath));
                lengths.push_back(settings.nodeWidth - node_event);
            }

            if (subtreeTime > parent.time)
            {
                branches.push_back(BranchPoint(&parent,this,0.0,false,onPath));
                lengths.push_back(min(time,subtreeTime)-parent.time);
                assert( lengths.back() > 0.0 );
            }

            return onPath;
        }
    }

    void GetAllLeaves(list<TreeNode*> &allLeaves, bool isRoot = true)
    {
        if (isLeaf)
            allLeaves.push_back(this);
        for (list<TreeNode*>::iterator i=children.begin(); i != children.end(); ++i)
        {
            (*i)->GetAllLeaves(allLeaves,false);
        }
    }


    // Choose a leaf at random to detach AND sample the Poisson events R,S,T (EXPERIMENTAL)
    TreeNode* GetRandomLeafAndSampleRST(BaseSettings &settings, double parentTime)
    {
        // sample R
        events.clear(); // delete old events

        if (isLeaf && settings.PoissonEventsOnFinalEdge>0) // if we are using the unbounded div. function method
        {
            for (int i=0; i<settings.PoissonEventsOnFinalEdge; i++) // sample L events from the prior on this edge
            {
                double At = settings.A(parentTime)+exp(settings.logDivergenceRateFactor(1))*(-log(settings.fRand(0,1)));
                double t= settings.invA(At);
                if (t<=parentTime || t>=time)
                    throw 1;
                events.push_back(t);
            }
        }
        else
        {
            double Aa = settings.Astar(parentTime) - settings.A(parentTime);
            double Ab = settings.Astar(time)- settings.A(time);
            double rate = (Ab-Aa)/exp(settings.logDivergenceRateFactor(numDescendants));
            int numEvents = settings.rpoiss(rate);
            cout << "numEvents: " << numEvents << endl;

            for (int i=0; i<numEvents; i++)
            {
                double t= settings.invAplusAstar(Aa+settings.fRand(0.0,1.0)*(Ab-Aa));
                if (t<=parentTime || t>=time)
                    throw 1;
                events.push_back(t);
            }
        }
        events.sort();

        if (isLeaf)
            return this;

        int temp=settings.iRand(0, numDescendants);
        int counter = 0;
        TreeNode* result = NULL;
        TreeNode* other=NULL;
        bool makeS = false;
        for (list<TreeNode*>::iterator i=children.begin(); i != children.end(); ++i)
        {
            counter+=(*i)->numDescendants;
            if (counter>=temp && (result == NULL))
            {
                if ((*i)->isLeaf && ListHasOnlyNElements<TreeNode*>(children,2))
                    makeS = true;
                result=(*i)->GetRandomLeafAndSampleRST(settings, time);
            }
            else
            {
                other = *i;
                (*i)->SampleT(settings,time);
            }
        }
        if (makeS) // if detaching will remove this node, it becomes S
        {
            if (other->isLeaf)
            {
                other->events.clear();
                for (int i=0; i<settings.PoissonEventsOnFinalEdge-1; i++) // sample L-1 events from the prior on this edge (since one will be S)
                {
                    double At = settings.A(parentTime)+exp(settings.logDivergenceRateFactor(1))*(-log(settings.fRand(0,1)));
                    double t= settings.invA(At);
                    if (t<=parentTime || t>=1.0)
                        throw 1;
                    other->events.push_back(t);
                }
            }
            else
                other->events.insert(other->events.begin(),events.begin(),events.end());
            // TODO why don't we have to add our divergence events here?
            other->events.push_back(time);
            other->events.sort();
        }
        return result;
    }


    // Choose a leaf at random to detach AND sample the Poisson events R,S,T (EXPERIMENTAL)
    // TODO: I don't think this works for unbounded divergence functions yet
    TreeNode* GetRandomSubtreeAndSampleRST(BaseSettings &settings, double parentTime, bool isRoot=true)
    {

        events.clear();
        TreeNode* result = NULL;
        // if we got to a leaf then we must choose this
        if (isLeaf)
            return this;
        // return this subtree?
        double probMe = settings.subtreeInflationFactor / (settings.subtreeInflationFactor + (double)numDescendants);
        if ((settings.fRand() < probMe) && !isRoot) // the subtree rooted at this node will be detached
            return this;

        // if not, randomly select a branch to go down
        int temp=settings.iRand(0, numDescendants-1);
        int counter = 0;
        // choose a branch at random, proportional to the number of descendants for each
        TreeNode* chosenBranch = NULL;
        for (list<TreeNode*>::iterator i=children.begin(); i != children.end(); ++i)
        {
            counter+=(*i)->numDescendants;
            if (counter>=temp) // choose this branch
            {
                chosenBranch=*i;
                break;
            }
        }
        result=chosenBranch->GetRandomSubtreeAndSampleRST(settings, time, false); // recurse down this branch

        // sample T on the remainder of the tree
        TreeNode* other=NULL;
        for (list<TreeNode*>::iterator i=children.begin(); i != children.end(); ++i)
        {
            if (chosenBranch != *i)
            {
                other = *i;
                (*i)->SampleT(settings,time,result->numDescendants,result->time);
            }
        }

        // sample R
        double Aa = settings.Astar(parentTime) - settings.A(parentTime);
        double Ab = settings.Astar(time)- settings.A(time);
        double rate = (Ab-Aa)/exp(settings.logDivergenceRateFactor(numDescendants-result->numDescendants,result->numDescendants));
        int numEvents = settings.rpoiss(rate);
        for (int i=0; i<numEvents; i++)
        {
            double t= settings.invAplusAstar(Aa+settings.fRand(0.0,1.0)*(Ab-Aa));
            if (t<=parentTime || t>=time)
                throw 1;
            if (t<result->time)
                events.push_back(t);
        }
        events.sort();

        // if this child is the result and we only have two children we will be deleted
        // when the result is removed, so our position becomes S _on our other child_
        // we also have to copy our events onto our other child since these will be lost otherwise
        if ((result==chosenBranch) && ListHasOnlyNElements<TreeNode*>(children,2)) // if detaching will remove this node, it becomes S _on our other child_
        {
            other->events.insert(other->events.begin(),events.begin(),events.end());
            other->events.push_back(time);
            other->events.sort();
        }

        return result;
    }

    // generate a list of events on the branch from parent to this (EXPERIMENTAL)
    void SampleT(BaseSettings &settings, double parentTime, int numInSubtree = 1, double subtreeTime= 1.0)
    {
        events.clear();

        if (parentTime < subtreeTime)
        {
            if (isLeaf && settings.PoissonEventsOnFinalEdge>0)
            {
                for (int i=0; i<settings.PoissonEventsOnFinalEdge; i++)
                {
                    double At = settings.A(parentTime)+exp(settings.logDivergenceRateFactor(1,numInSubtree))*(-log(settings.fRand()));
                    double t= settings.invA(At);
                    if (t<=parentTime || t>=time)
                        throw 1;
                    if (t<subtreeTime)
                        events.push_back(t);
                }
            }
            else
            {
                double Aa = settings.Astar(parentTime);
                double Ab = settings.Astar(time);
                double rate = (Ab-Aa)/exp(settings.logDivergenceRateFactor(numDescendants,numInSubtree));
                int numEvents = settings.rpoiss(rate);
                for (int i=0; i<numEvents; i++)
                {
                    double t= settings.invAstar(Aa+settings.fRand()*(Ab-Aa));
                    if (t<=parentTime || t>=time)
                        throw 1;
                    if (t<subtreeTime)
                        events.push_back(t);
                }
            }
        }

        events.sort();

        for (list<TreeNode*>::iterator i=children.begin(); i != children.end(); ++i)
        {
            (*i)->SampleT(settings, time, numInSubtree, subtreeTime);
        }
    }

    void ClearEvents()
    {
        events.clear();
        for (list<TreeNode*>::iterator i=children.begin(); i != children.end(); ++i)
        {
            (*i)->ClearEvents();
        }
        node_event=0.0;
    }

    double LatestInternalNode()
    {
        if (isLeaf)
            return 0.0;
        else
        {
            double m=time;
            for (list<TreeNode*>::iterator i=children.begin(); i != children.end(); ++i)
            {
                m=max((*i)->LatestInternalNode(),m);
            }
            return m;
        }
    }

    // TODO: make this multithreaded. You would need a little class which
    // stores the input argument and the resulting branchPoints and logProbs. Once executed you
    // could collect up all the branchPoints and logProbs of your children. 
    void EvaluatePotentialBranchPoints(std::vector<BranchPoint> &branchPoints,
                                       list<double> &logProbs,
                                       TreeNode* subtree,
                                       TreeNode &parent,
                                       double logProbOfGettingHere,
                                       double logProdOneMinusD,
                                       BaseSettings &s)
    {
        double checksorted = 0.0;
        for (list<double>::iterator i=events.begin(); i != events.end(); ++i)
        {
            if (*i < checksorted)
                throw 1;
            checksorted=*i;
            BranchPoint bp(&parent,this,*i,false);
            branchPoints.push_back(bp);
            // likelihood calculation
            GaussianVector msg[3];
            double temp = *i;

            msg[0]=SampleAverageConditional(subtree->Marg_Location, subtree->time - *i);
            msg[1]=SampleAverageConditional(Marg_Location, time - *i);
            msg[2]=SampleAverageConditional(parent.Marg_Location, *i - parent.time);
            GaussianVector prod = msg[0] * msg[1] * msg[2];
            double pY = prod.GetLogNormalizer();
            for (int j=0; j<3; j++)
                pY -= msg[j].GetLogNormalizer();
            assert(!isnan<double>(pY));

            // p(x_b|x_a) normalisation constant
            GaussianVector tob = SampleAverageConditional(parent.Marg_Location, time - parent.time);
            pY -= sumVector(tob.GetLogProb(Marg_Location.GetMean()));
            if (isnan<double>(pY))
                throw 1;
            // no structure contribution, already accounted for in Poisson process
            // probability of s-1 datapoints following this branch point
            if (subtree->numDescendants>1)
            {
                pY += lgamma((double)subtree->numDescendants-s.alpha) - lgamma(1.0-s.alpha);
                pY -= lgamma((double)subtree->numDescendants+(double)numDescendants+s.theta) - lgamma((double)numDescendants+s.theta+1.0);
            }
            assert(!isnan<double>(pY));
            // include probability of branching to here
            pY += logProbOfGettingHere;
            // if appropriate include contribution from d(i)
            if (!(isLeaf && s.PoissonEventsOnFinalEdge>0))
            {
                // should there be a ratio of gamma terms here for diverging?
                double logdu = log(s.a(*i))-log(s.astar(*i));
                // is it correct that is constant? probably...
                //cout << logdu << endl;
                pY += logdu;
                pY += logProdOneMinusD;
                logProdOneMinusD += log(1.0-exp(logdu));
            }
            else
            {
                // pY += log(s.a(*i)); // should there be a ratio of gamma terms here for diverging?
                pY += logProdOneMinusD;
                pY -= log((double)s.PoissonEventsOnFinalEdge);
            }
            assert(!isnan<double>(pY));
            logProbs.push_back(pY);
            //cout << "Evaluating branch at t=" << temp << " pY=" << pY << endl;
        }

        if (isLeaf)
            return;
        if (time>=subtree->time)
            return;
        // add possibility of forming a new branch at this existing branch point
        BranchPoint ebp(&parent, this, time, true);
        branchPoints.push_back(ebp);
        // simple likelihood calculation
        GaussianVector down = SampleAverageConditional(Marg_Location, subtree->time - time);
        double probY = sumVector(down.GetLogProb(subtree->Marg_Location.GetMean())); // TODO should we try to integrate over x_u here?
        double probNewBranch = log(s.theta+s.alpha*(double)children.size()) - log((double)numDescendants+s.theta);
        // acount for the fact that the remaining s-1 datapoints are forced to go down this branch
        if (subtree->numDescendants>1)
        {
            probNewBranch += lgamma((double)subtree->numDescendants-s.alpha) - lgamma(1.0-s.alpha);
            probNewBranch -= lgamma((double)subtree->numDescendants+(double)numDescendants+s.theta) - lgamma((double)numDescendants+s.theta+1.0);
        }
        // TODO: shouldn't Poisson process deal with these events? Need thinning?
        double pY = logProbOfGettingHere + probY + probNewBranch;
        cout << "Log prob new branch: " << pY << endl;
        //logProbOfGettingHere -= log(1.0-exp(probNewBranch));
        if (isnan<double>(pY))
            throw 1;
        if (exp(pY)>0.0)
            logProbs.push_back(pY);
        //boost::thread_group tg;
        for (list<TreeNode*>::iterator i=children.begin(); i != children.end(); ++i)
        {
            double logProbThisBranch = logProbOfGettingHere;
            if (subtree->numDescendants==1) // probability of 1 datapoint going down this branch
            {
                logProbThisBranch += (log((double)(*i)->numDescendants-s.alpha) - log((double)numDescendants+s.theta));
            }
            else // probability of s datapoints going down this branch
            {
                // NB could calculate these without lgamma but probably (?) cheaper this way, certainly for larger subtrees
                logProbThisBranch += (lgamma((double)(*i)->numDescendants+(double)subtree->numDescendants-s.alpha)-lgamma((double)(*i)->numDescendants-s.alpha));
                logProbThisBranch -= (lgamma((double)numDescendants+(double)subtree->numDescendants+s.theta)-lgamma((double)numDescendants+s.theta));
            }
            (*i)->EvaluatePotentialBranchPoints(branchPoints, logProbs, subtree, *this, logProbThisBranch, logProdOneMinusD, s);
            //tg.create_thread(&TreeNode::EvaluatePotentialBranchPoints, *i, branchPoints, logProbs, subtree, *this, logProbThisBranch, logProdOneMinusD, s);
        }
        //tg.join_all();
    }
};

// attach subtree at this branch point and maybe resample parent location
void BranchPoint::AttachSubtree(TreeNode* subtree, BaseSettings& settings)
{
    if (isExistingBranchPoint) // this is a new branch at an existing branch point
    {
        child->numDescendants+=subtree->numDescendants;
        child->children.push_back(subtree);
	subtree->parent=child; 
	subtree->Msg_Normal_ParentLocation = SampleAverageConditional(subtree->Marg_Location / subtree->Msg_Normal_Location, subtree->time - child->time); 
	child->Marg_Location *= subtree->Msg_Normal_ParentLocation; 
	child->bpSweepUp(subtree->numDescendants);
	
	subtree->bpSweepDown(*child); // this will update subtree->Msg_Normal_Location
        // optionally we resample the location of the parent since the addition of the subtree may have a
        // large effect on its location
        if (settings.ResampleAttachmentLocation){
            GaussianVector marg=SampleAverageConditional(parent->Marg_Location, child->time - parent->time);
            for (list<TreeNode*>::iterator i=child->children.begin(); i != child->children.end(); ++i)
	      marg *= SampleAverageConditional((*i)->Marg_Location, (*i)->time-child->time);
	    child->Marg_Location.SetPoint(marg.Sample(settings.gen));
        }
    }
    else
    {
        TreeNode *new_node = new TreeNode(time,false);
        new_node->children.push_back(child);
	child->parent=new_node; 
        new_node->children.push_back(subtree);
	subtree->parent=new_node; 
        new_node->numDescendants = child->numDescendants+subtree->numDescendants;
        parent->children.remove(child);
        parent->children.push_back(new_node);
	new_node->parent=parent; 
        
        // to continue sampling subtree locations we need to sample a location for this new node
	// TODO: don't do this if integrating! 
	if (0){
	  parent->numDescendants+=subtree->numDescendants;
	  GaussianVector msg[3];
	  msg[0]=SampleAverageConditional(subtree->Marg_Location, subtree->time - time);
	  msg[1]=SampleAverageConditional(child->Marg_Location, child->time - time);
	  msg[2]=SampleAverageConditional(parent->Marg_Location, time - parent->time);
	  GaussianVector prod = msg[0] * msg[1] * msg[2];
	  new_node->Marg_Location.SetPoint(prod.Sample(settings.gen));
	} else {
	  // oh god there are a lot of connections to get right here...
	  GaussianVector wp_to_norm = parent->Marg_Location / child->Msg_Normal_ParentLocation; 
	  new_node->Msg_Normal_Location = SampleAverageConditional(wp_to_norm, new_node->time - parent->time);
	  GaussianVector wc_to_norm = child->Marg_Location / child->Msg_Normal_Location; 
	  subtree->Msg_Normal_ParentLocation = SampleAverageConditional(subtree->Marg_Location / subtree->Msg_Normal_Location, subtree->time - new_node->time);
	  child->Msg_Normal_ParentLocation = SampleAverageConditional(wc_to_norm, child->time - new_node->time); 
	  new_node->Marg_Location = new_node->Msg_Normal_Location * subtree->Msg_Normal_ParentLocation * child->Msg_Normal_ParentLocation;
	  new_node->Msg_Normal_ParentLocation = SampleAverageConditional(subtree->Msg_Normal_ParentLocation * child->Msg_Normal_ParentLocation, new_node->time - parent->time); 
	  child->Marg_Location /= child->Msg_Normal_Location; 
	  child->Msg_Normal_Location = SampleAverageConditional(new_node->Msg_Normal_Location * subtree->Msg_Normal_ParentLocation, child->time - new_node->time);
	  child->Marg_Location *= child->Msg_Normal_Location;
	  subtree->bpSweepDown(*new_node); 
	  parent->numDescendants += subtree->numDescendants; 
	  parent->bpSweepUp(subtree->numDescendants);
	}
    }
};

string BranchPoint::ToString()
{

    return "Time " +  boost::lexical_cast<string>(time) + " parent " + boost::lexical_cast<string>(parent->time) + " child " + boost::lexical_cast<string>(child->time) + " onPath " + boost::lexical_cast<string>(onPath) + " isExistingBP " + boost::lexical_cast<string>(isExistingBranchPoint);

}

bool BranchPoint::IsBefore(BranchPoint &other)
{
    if (isExistingBranchPoint && other.isExistingBranchPoint)
    {
        if (child->time == other.child->time)
            return time < other.time;
        else
            return child->time < other.child->time;
    }
    else
    {
        double myRealTime = isExistingBranchPoint ? child->time : time;
        double otherRealTime = other.isExistingBranchPoint ? other.child->time : other.time;
        return myRealTime < otherRealTime;
    }
}

#endif
