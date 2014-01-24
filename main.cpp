#define _USE_MATH_DEFINES


#include <iostream>
#include <fstream>
#include <stdio.h>
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
//#include "TreeNode.hpp"
#include "kingman.hpp"
#include <string>
#include <limits>
#include <algorithm>
#include <assert.h>

#include "Tree.hpp"
#include "pydtsamplers.hpp"

using namespace std;


GaussianVector stdToBlasHelper(std::vector<double> &x)
{
    GaussianVector retval(x.size());
    double prec=100;
    for (int i=0; i<x.size(); i++)
    {
        if (x[i]==0.0)
        {
            retval.MeanTimesPrecision[i]=0.0;
            retval.Precision[i]=0.0;
        }
        else
        {
            retval.MeanTimesPrecision[i] = prec * x[i];
            retval.Precision[i] = prec;
        }
    }
    return retval;
}

boost::numeric::ublas::vector<double> stdToBlas(std::vector<double> &x)
{
    boost::numeric::ublas::vector<double> retval(x.size());

    for (int i=0; i<x.size(); i++)
        retval[i]=x[i];
    return retval;
}


std::vector<double> readRow(std::string row) {
    std::vector<double> retval;
    std::istringstream is(row);
    double num;
    while (is >> num)
        retval.push_back(num);

    return retval;
}

std::vector<boost::numeric::ublas::vector<double> > readVector(std::istream &is) {
    std::string line;
    std::vector<boost::numeric::ublas::vector<double> > retval;
    while (std::getline(is, line))
    {
        std::vector<double> b=readRow(line);
        boost::numeric::ublas::vector<double> a;
        a=stdToBlas(b);
        retval.push_back(a);
    }
    return retval;
}

std::vector<string > readNames(char *filename)
{
    std::vector<string > retval;
    ifstream in(filename);
    std::string line;
    string st;
    while (std::getline(in, line))
    {
        std::istringstream is(line);
        is >> st;
        retval.push_back(st);
    }
    in.close();
    return retval;
}

std::vector<boost::numeric::ublas::vector<double> > readFile(char *filename)
{
    ifstream in(filename);
    std::vector<boost::numeric::ublas::vector<double> > retval = readVector(in);
    in.close();
    return retval;
}


void make_train_test(char* datafile,  int seed)
{
    // read data
    std::vector<boost::numeric::ublas::vector<double> > data=readFile(datafile);

    int totalN = data.size();
    int testN = totalN / 10;
    int trainN = totalN - testN;
    std::list<boost::numeric::ublas::vector<double> > test_data(testN);
    std::vector<boost::numeric::ublas::vector<double> > train_data(trainN);

    std::random_shuffle ( data.begin(), data.end() ); // should also shuffle the labels!

    std::copy(data.begin(),data.begin()+trainN,train_data.begin());
    std::copy(data.begin()+trainN,data.end(),test_data.begin());

    ofstream outputfile;
    string fn = "splits/train" + boost::lexical_cast<string>(seed);
    outputfile.open(fn.c_str(),ios::out);
    for (std::vector<boost::numeric::ublas::vector<double> >::iterator i=train_data.begin(); i != train_data.end(); ++i)
        outputfile << (*i) << endl;
    outputfile.close();

    string fn2 = "splits/test" + boost::lexical_cast<string>(seed);
    outputfile.open(fn2.c_str(),ios::out);
    for (std::list<boost::numeric::ublas::vector<double> >::iterator i=test_data.begin(); i != test_data.end(); ++i)
        outputfile << (*i) << endl;
    outputfile.close();

}

void data_test(bool mh, double theta, char* datafile, char* labelsFile, int iterations, bool flat, int seed)
{
    // read data
    std::vector<boost::numeric::ublas::vector<double> > data=readFile(datafile);

    int totalN = data.size();
    int testN = totalN / 10;
    int trainN = totalN - testN;
    std::list<boost::numeric::ublas::vector<double> > test_data(testN);
    std::vector<boost::numeric::ublas::vector<double> > train_data(trainN);

    std::random_shuffle ( data.begin(), data.end() ); // should also shuffle the labels!

    std::copy(data.begin(),data.begin()+trainN,train_data.begin());
    std::copy(data.begin()+trainN,data.end(),test_data.begin());

    // read names of objects
    std::vector<string > names = readNames(labelsFile);
    int numLeaves = data.size();

    Settings settings;
    settings.D=data[0].size();

    if (theta != -1.0)
    {
        settings.sampleTheta=false;
        settings.sampleAlpha=false;
        settings.theta=theta;
    }

    //TreeNode* root = new TreeNode(0.1,false);
    TreeNode* firstLeaf = new TreeNode(1,true);
    firstLeaf->label=names[0];
    //firstLeaf->Marg_Location.resize(settings.D);
    firstLeaf->Marg_Location.SetPoint(train_data[0]);
    //root->children.push_back(firstLeaf);
    list<TreeNode*> leaves;
    leaves.push_back(firstLeaf);
    for (int i=1; i<trainN; i++)
    {
        TreeNode* leaf = new TreeNode(1.0, true);
        leaf->label=names[i];
        //leaf->Marg_Location.resize(settings.D);
        leaf->Marg_Location.SetPoint(train_data[i]);
        leaves.push_back(leaf);
    }

    Settings binary;
    binary.D= settings.D;
    binary.theta=0.0;
    binary.alpha=0.0;

    Tree*  learntTree =  flat ? Tree::FlatTree(leaves,settings) : Tree::RandomTree(leaves,binary);

    cout << "initial tree has " << learntTree->countLeaves(true) << " leaves and " << learntTree->numInternalNodes() << " internal nodes " << endl;


    if (mh)
    {
        std::string fn = std::string("cluster_results/mh_") + boost::lexical_cast<string>(theta) + (flat ? "_flat" : "_rand" ) + boost::lexical_cast<string>(seed);
        learntTree = MH(learntTree, test_data, settings, iterations, (fn + "_ml.txt").c_str(), (fn + "_newick.txt").c_str());
    }
    else
    {
        std::string fn = std::string("cluster_results/slice_") + boost::lexical_cast<string>(theta)+ (flat ? "_flat" : "_rand" ) + boost::lexical_cast<string>(seed); ;
        learntTree = SliceSample(learntTree, test_data, settings, iterations, (fn + "_ml.txt").c_str(), (fn + "_newick.txt").c_str());
    }

    cout << "ML learnt tree: " << learntTree->MarginalLikelihood(settings) << endl;


    delete learntTree;
}

void test_sampler(bool mh)
{
    Settings settings;
    settings.D=2;
    settings.theta = 0.0;
    TreeNode* root = new TreeNode(0.5,false);
    TreeNode* firstLeaf = new TreeNode(1,true);
    firstLeaf->label="0";
    root->children.push_back(firstLeaf);
    TreeNode* secondLeaf = new TreeNode(1,true);
    secondLeaf->label="1";
    root->children.push_back(secondLeaf); 
    list<TreeNode*> leaves;
    Tree* tree=new Tree(root, settings);
    leaves.push_back(firstLeaf);
    leaves.push_back(secondLeaf);
    int numLeaves = 10; 
    for (int i=0; i<numLeaves-2; i++)
    {
        TreeNode* leaf = new TreeNode(1.0, true);
        leaf->label=boost::lexical_cast<string>(i+1);
        tree->AddChild(settings, leaf);
        leaves.push_back(leaf);
    }

    tree->sampleData(settings);
    cout << tree->newick() << endl;
    cout << "ML true tree: " << tree->MarginalLikelihood(settings) << endl;

    //return;

    //Tree* flatTree = Tree::FlatTree(leaves,settings);

    //BoundedSettings boundedSettings;
    //boundedSettings.D=settings.D;

    //cout << "ML true tree (bounded a(t)): " << tree->MarginalLikelihood(boundedSettings) << endl;

    //cout << "ML flat tree: " << flatTree->MarginalLikelihood(settings) << endl;

    //Tree* builtTree = BuildTree(leaves, settings);

    //cout << "ML built tree: " << builtTree->MarginalLikelihood(settings) << endl;
    // return ;

    Tree*  learntTree = Tree::RandomTree(leaves,settings);

    cout << "initial tree has " << learntTree->countLeaves(true) << " leaves" << endl;

    list<boost::numeric::ublas::vector<double> > dummy;
    //learntTree = PoissonProcessSamplerWrapper(learntTree, settings, 100, 10, "ps_results.txt", dummy);

    if (mh)
        learntTree = MH(learntTree, dummy, settings, 100000, "mh_results.txt");
    else
        learntTree = SliceSample(learntTree, dummy, settings, 10000, "slice_results.txt");

    cout << "ML learnt tree: " << learntTree->MarginalLikelihood(settings) << endl;

    delete tree;
    delete learntTree;
}


void gettingItRightMH(bool posterior, int repeats = 1000, int nleaves = 10)
{
    Settings settings;
    settings.D=1;
    //settings.alpha=-.5;
    //settings.theta=-2.0*settings.alpha;
    settings.theta=1.0;
    TreeNode* firstLeaf = new TreeNode(1.0,true);
    firstLeaf->label="0";
    list<TreeNode*> leaves;
    Tree* tree=new Tree(firstLeaf, settings);
    leaves.push_back(firstLeaf);
    settings.c=1;
    for (int i=1; i<nleaves; i++)
    {
        TreeNode* leaf = new TreeNode(1.0, true);
        leaf->label=boost::lexical_cast<string>(i);
        tree->AddChild(settings, leaf);
        //cout << tree->countLeaves() << endl;
        leaves.push_back(leaf);
    }
    //cout << tree->newick() << endl;
    tree->sampleData(settings);
    settings.c=1;
    //cout << tree->newick() << endl;
    char *fn = "ps_results.txt";
    //remove(fn);
    list<boost::numeric::ublas::vector<double> > dummy;
    if (posterior)
    {
        for (int i=1; i<repeats; i++)
        {
            tree->sampleData(settings);

            //cout <<  (*tree->root->children.begin())->Marg_Location.GetMean()[0] - (*++tree->root->children.begin())->Marg_Location.GetMean()[0] << endl;
            //tree = MH(tree, settings, 10, "mh_results.txt");
            //list<boost::numeric::ublas::vector<double> > test;

            tree = SliceSample(tree, dummy, settings, 10, fn);
            cout <<  tree->root->time << endl;
            //cout << tree->root->countLeaves(false, true) << endl;
        }
    }
    else
    {
        for (int i=1; i<repeats; i++)
        {
            Tree* randomTree = Tree::RandomTree(leaves,settings);
            //randomTree->sampleData(settings);
            //cout <<  (*randomTree->root->children.begin())->Marg_Location.GetMean()[0] - (*++randomTree->root->children.begin())->Marg_Location.GetMean()[0] << endl;
            cout <<  randomTree->root->time << endl;
            //cout << randomTree->root->countLeaves(false, true) << endl;
            delete randomTree;
        }
    }
//
    //cout << "ML random tree: " << randomTree->MarginalLikelihood(settings) << endl;


//   Tree*  learntTree = randomTree;
//   cout << "initial tree has " << learntTree->countLeaves(true) << " leaves" << endl;
//   list<boost::numeric::ublas::vector<double> > test;
//   //learntTree = PoissonProcessSamplerWrapper(learntTree, settings, 10000, 100, "ps_results.txt", test);
//   learntTree = MH(learntTree, settings, 1, "mh_results.txt");
//
//   cout << "ML learnt tree: " << learntTree->MarginalLikelihood(settings) << endl;
//
    delete tree;
    //delete randomTree;
    //delete learntTree;
}

void dataTest(list<TreeNode*> leaves, list<boost::numeric::ublas::vector<double> > &test, int D)
{
    Settings settings;
    settings.D=D;


    Tree* randomTree = Tree::RandomTree(leaves,settings);
    cout << "ML random tree: " << randomTree->MarginalLikelihood(settings) << endl;
    //BoundedSettings boundedSettings;
    //boundedSettings.D=settings.D;

    //cout << "ML true tree (bounded a(t)): " << tree->MarginalLikelihood(boundedSettings) << endl;
    //Tree* flatTree = Tree::FlatTree(leaves,settings);
    //cout << "ML flat tree: " << flatTree->MarginalLikelihood(settings) << endl;

    Tree* builtTree = BuildTree(leaves, settings);

    cout << "ML built tree: " << builtTree->MarginalLikelihood(settings) << endl;
    // return ;

    Tree*  learntTree = randomTree;

    cout << "initial tree has " << learntTree->countLeaves(true) << " leaves" << endl;

    learntTree = PoissonProcessSamplerWrapper(learntTree, settings, 100, 100, "ps_results.txt", test);
    //learntTree = MH(learntTree, settings, 1000, "mh_results.txt");
    learntTree->OutputTree("tree.txt");
    cout << "ML learnt tree: " << learntTree->MarginalLikelihood(settings) << endl;

    // delete learntTree;
    delete builtTree;
    delete randomTree;
}

// Sample "repeats" trees with "treeSize" leaves, and output a histogram of the number
// of leaves at each depth
void tree_hist(Settings &settings, int treeSize = 1000, int repeats = 1000)
{
    map<int,int> hist;

    for (int rep=0; rep<repeats; rep++)
    {
        TreeNode root(0.1,false);
        TreeNode *firstLeaf = new TreeNode(1,true);
        firstLeaf->label="0";
        root.children.push_back(firstLeaf);
        for (int i=0; i<treeSize; i++)
        {
            TreeNode* leaf = new TreeNode(1.0, true);
            leaf->label=boost::lexical_cast<string>(i+1);
            root.AddChild(settings, leaf);
            //cout << root.countLeaves() << endl;
        }
        root.LeafDepthHist(hist);
    }
    for (int i=0; i<300; i++)
        cout << i << " " << (hist.count(i)>0 ? hist[i] : 0) << endl;
}


// Sample "repeats" trees with "treeSize" leaves, and output a histogram of the number
// of leaves at each depth
void tree_hist_kingmans(Settings &settings, int treeSize = 1000, int repeats = 1000)
{
    map<int,int> hist;

    for (int rep=0; rep<repeats; rep++)
    {
        TreeNode* root = KingmansCoalescent(treeSize, settings);
        root->LeafDepthHist(hist);
        delete root;
    }
    for (int i=0; i<300; i++)
        cout << i << " " << (hist.count(i)>0 ? hist[i] : 0) << endl;
}

void tree_hist_single(Settings &settings, int treeSize = 1000, int repeats = 100000)
{
    map<int,int> hist;

    for (int rep=0; rep<repeats; rep++)
    {
        TreeNode root(0.1,false);
        TreeNode *firstLeaf = new TreeNode(1,true);
        firstLeaf->label="0";
        root.children.push_back(firstLeaf);
        for (int i=0; i<treeSize; i++)
        {
            TreeNode* leaf = new TreeNode(1.0, true);
            leaf->label=boost::lexical_cast<string>(i+1);
            root.AddChild(settings, leaf);
            //cout << root.countLeaves() << endl;
        }
        map<int,int> treehist;
        root.LeafDepthHist(treehist);
        // otal = std::accumulate(treehist.begin(), treehist.end(), 0)
        int h = settings.sampleMap(treehist,treeSize);
        if (hist.count(h)==0)
            hist[h]=0;
        hist[h]+=1;
    }
    for (int i=0; i<300; i++)
        cout << i << " " << (hist.count(i)>0 ? hist[i] : 0) << endl;
}

void depth_as_func_of_n(Settings &settings, int maxN = 1e4, int repeats=100)
{
    double depths [maxN];
    fill(depths,depths+maxN,0);
    for (int r=0; r<repeats; r++)
    {
        TreeNode root(0.1,false);
        TreeNode *firstLeaf = new TreeNode(1,true);
        firstLeaf->label="0";
        root.children.push_back(firstLeaf);
        for (int i=0; i<maxN; i++)
        {
            TreeNode* leaf = new TreeNode(1.0, true);
            leaf->label=boost::lexical_cast<string>(i+1);
            depths[i]+=(double)root.AddChild(settings, leaf)/(double)repeats;
        }
    }

    for (int i=0; i<maxN; i++)
        cout << i << " " << depths[i] << endl;
}

int main (int argc,char *argv[]) {

    //Settings settings;
     // if (argc < 6)
     // {
     //     cout << "Not enough inputs, usage is" << endl;
     //     cout << "ts mh(rather than slice) theta datafile namesfile iterations flat(rather than random) seed" << endl;
     //     return 1;
     // }
     
     //int seed = boost::lexical_cast<int>(argv[7]);
   // std::srand ( seed );
  /*    for (int seed=1; seed <= 10; seed++)
    {
        std::srand ( seed );
        make_train_test(argv[1], seed);
	}*/
// gettingItRightMH(boost::lexical_cast<bool>(argv[1]), boost::lexical_cast<int>(argv[2]), boost::lexical_cast<int>(argv[3]));
    //settings.theta=boost::lexical_cast<double>(argv[1]);
    //settings.alpha=boost::lexical_cast<double>(argv[2]);
    //settings.c=boost::lexical_cast<double>(argv[3]); */

    //depth_as_func_of_n(settings);
    //testMH();
    test_sampler(boost::lexical_cast<bool>(argv[1]));
    // data_test(boost::lexical_cast<bool>(argv[1]),boost::lexical_cast<double>(argv[2]),argv[3],argv[4],boost::lexical_cast<int>(argv[5]),boost::lexical_cast<bool>(argv[6]), seed);
    //TreeNode *root = KingmansCoalescent(10, settings);
    //tree_hist_kingmans(settings);

    return 0;
}
