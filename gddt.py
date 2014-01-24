# -*- coding: utf-8 -*-
import sys
from numpy import *
import time
import pylab
import logging
from scipy.special import gammaln, gamma
from scipy.stats import rv_discrete
import pydot
G=pydot.Dot()
NODECOUNTER=1

def randsample(p):
    p=rv_discrete(name='adhoc', values=(range(len(p)), p))
    return p.rvs()
    
def H(n,a,b):
    i=arange(1,n+1)
    return sum( exp( gammaln(i-b) - gammaln(i+1+a) ))
    
def J(n,a,b):
    res=H(sum(n),a,b)
    for i in n:
        res-=H(i,a,b)
    return res

PLOT_INTERNAL_NODES=False
PLOT_TREE=False

alpha=0.0
beta=0.0
C=1.0
col="gbrmcyk"
class node:
    __slots__=["children","location","time","At","count"]
    
    def __init__(self,time,At):
        self.children=[]
        self.location=None
        self.time=time
        self.count=[]
        self.color=1
        self.At=At
        
    def isLeaf(self):
        return self.children==[]

    def newick(self,names,parenttime=0.0):
        if self.isLeaf():
            #print self.label
            return "%s:%f"%(names[self.label-1],self.time-parenttime)
        else:
            res="("
            for i in range(len(self.children)):
                res+=self.children[i].newick(names,self.time)
                if i<(len(self.children)-1):
                    res+=","
            return res+"):%f"%self.time
    
    def addPoint(self,invA,branchLabel=1):
        if len(self.children)>1:
            prob=zeros(len(self.count)+1)
            prob[0:len(self.count)]=array(self.count)-beta
            prob[len(self.count)]=alpha+beta*len(self.count)
            prob=prob/prob.sum()
            i=randsample(prob)
        else:
            i=0
        if i==len(self.count): # new block
            #print "New blocK!"
            leaf=node(1,inf)
            leaf.color=branchLabel
            self.children.append(leaf)
            self.count.append(1)
        else: # block i
            At=self.At+exp(gammaln(self.count[i]+1+alpha)-gammaln(self.count[i]-beta))*-log(random.rand())
            #print At
            td=invA(At)
            #print td
            if td>self.children[i].time:
                self.children[i].addPoint(invA,branchLabel)
            else:
                temp=node(td,At)
                temp.color=self.children[i].color
                temp.children.append(self.children[i])
                temp.count.append(self.count[i])
                leaf=node(1,inf)
                leaf.color=branchLabel
                temp.children.append(leaf)
                temp.count.append(1)
                self.children[i]=temp
            self.count[i]+=1
            assert td<=1.0, "td=%f"%td
            
    def ncrpAddPoint(self,A,S):
        
        prob=zeros(len(self.count)+1)
        prob[0:len(self.count)]=array(self.count)-beta
        prob[len(self.count)]=A(self.time)+beta*len(self.count)
        #print A(self.time)
        prob=prob/prob.sum()
        i=randsample(prob)
        if i==len(self.count): # or self.children[i].isLeaf(): # new block
            #print "new block"
            self.ncrpDiverge(S)
            self.label.set_width(.5)
        else: # block i
            #print "old block"
            self.children[i].ncrpAddPoint(A,S)
            self.count[i]+=1
        
    def diffuse(self,parentLocation,parentTime,sigma):
        self.location=parentLocation+sigma*(self.time-parentTime)*random.randn(len(parentLocation))
        for c in self.children:
            c.diffuse(self.location,self.time,sigma)
            
    def diffusePlot(self,parentLocation,parentTime,sigma,delta=0.001):
        #assert len(parentLocation)==1, "Dimension should be 1"
        if self.time!=parentTime:
            t=arange(parentTime,self.time,delta)
            #print t
            x=random.randn(len(t))*delta
            x[0]=0.0
            x=parentLocation+x.cumsum()
            pylab.plot(t,x,col[self.color])
            print "Here"
            self.location=x[-1]
        else:
            self.location=parentLocation
        pylab.scatter(self.time,self.location)
        for c in self.children:
            c.diffusePlot(self.location,self.time,sigma)
            
    def plot(self,N=None,power=1.0):
        if N==None:
            N=self.countChildren()
        if self.isLeaf():
            pylab.scatter(self.location[0],self.location[1],c='k',alpha=.2)
        else:
            if PLOT_INTERNAL_NODES:
                pylab.scatter(self.location[0],self.location[1],c='g')
            if PLOT_TREE:
                for i in range(len(self.children)):
                    pylab.plot([self.location[0],self.children[i].location[0]],[self.location[1],self.children[i].location[1]],alpha=(float(self.count[i])/N)**power,c='k')
          
        for c in self.children:
            c.plot(N,power=power)
    
    def plotTime(self,N=None,gradient=True):
        if N==None:
            N=sum(self.count)
        alpha=1.0
        for i in range(len(self.children)):
            if gradient:
                alpha=float(self.count[i])/N
            pylab.plot([self.time,self.children[i].time],[self.location[0],self.children[i].location[0]],c='k',alpha=alpha)
            self.children[i].plotTime(N,gradient)
            
    def getDim(self,ind,points):
        if self.isLeaf():
            points.append(self.location[ind])
        for c in self.children:
            c.getDim(ind,points)
            
    def labelLeaves(self,counter=[1]):
        if self.isLeaf():
            #print counter
            self.label=counter[0]
            counter[0]+=1
        for c in self.children:
            c.labelLeaves(counter)
    
    def countChildren(self,checkConsistency=False):
        for i in range(len(self.children)):
            count=self.children[i].countChildren(checkConsistency)
            if checkConsistency:
                assert self.count[i]==count
            if len(self.count)<len(self.children):
                self.count=zeros(len(self.children))
            self.count[i]=count
        if self.isLeaf():
            return 1
        else:
            assert len(self.children)>=2
            return sum(self.count)
    
    def structureProbability(self):
        if self.isLeaf():  
            return 0
        K=len(self.count)
        if beta==0:
            result=sum(gammaln(self.count))-gammaln(sum(self.count)+alpha)
            if K>2:
                result+=(K-2)*log(alpha)
            if isnan(result):
                print self.count
                print alpha
        else:
            result=sum(gammaln(array(self.count)-beta))-gammaln(sum(self.count)+alpha)+(K-2)*log(beta)+gammaln(alpha/beta+K)-gammaln(alpha/beta+2)
        for c in self.children:
            result+=c.structureProbability()
        return result
        
    def timesProbability(self):
        if self.isLeaf():  
            return 0
        result = log(C) + (C * J(self.count,alpha,beta) - 1.0) * log(1.0-self.time)
        for c in self.children:
            result+=c.timesProbability()
        return result
    
    def getDataMatrix(self):
        dims=[]
        for i in range(0,len(self.location)):
            points=[]
            self.getDim(i,points)
            dims.append(points)
        return array(dims)
        
    def ncrpDiverge(self,S):
	if self.time>=1.0:
	  return
        newtime=self.time+1.0/double(S)
        newtime=min(1.0,newtime)
        newnode=node(newtime,0)
        global NODECOUNTER
        newnode.label=pydot.Node(name=str(NODECOUNTER),label=" ",width=.2,regular=True,fixedsize=True)
        NODECOUNTER=NODECOUNTER+1
        #G.add_edge(pydot.Edge(self.label,newnode.label))
        G.add_node(newnode.label)
        G.add_edge(pydot.Edge(self.label,newnode.label))
        #print self.label, newnode.label
        self.children.append(newnode)
        self.count.append(1)
        newnode.ncrpDiverge(S)
            
    def __str__(self):
        return ""
        
def DPMsample(N):
    clusters=[]
    counts=[]
    x=zeros([N,2])
    beta=0
    alpha=2
    for n in range(N):
        prob=zeros(len(clusters)+1)
        prob[0:len(clusters)]=array(counts)-beta
        prob[len(clusters)]=alpha+beta*len(clusters)
        prob=prob/prob.sum()
        i=randsample(prob)
        if i==len(clusters): 
            clusters.append(random.randn(2)*sqrt(3.0/4))
            counts.append(0)
        counts[i]=counts[i]+1
        x[n,:]=clusters[i]+random.randn(2)*sqrt(1.0/4)
    print counts
    pylab.scatter(x[:,0],x[:,1],color="k",alpha=.2)
    pylab.show()
    
def sampleDDT(N, D, sigma, invA):
    print "Creating tree"
    root=node(0,0)
    root.children.append(node(1,inf))
    root.count.append(1)
    for i in range(0,N-1):
        root.addPoint(invA,i+2)
    print "Count"
    root.children[0].countChildren(True)
    print "Prob"
    print root.children[0].structureProbability()
    print "Diffuse"
    root.diffuse(zeros(D),0,sigma)
    return root


def sampleNCRP(N, S, A):
    print "Creating tree"
    root=node(0,0)
    root.label=pydot.Node(name=0,label=" ",width=.1,regular=True)
    G.add_node(root.label)
    root.ncrpDiverge(S)
    for i in range(0,N-1):
        root.ncrpAddPoint(A,S)
    
    return root
    
def test(N=200,D=2):
    pylab.close('all')
    invA=lambda x: 1.0-exp(-x/C)
    print "Creating ddt"
    r=sampleDDT(N,D,1.0,invA)
    r.plot(N)
    #pylab.figure()
    print "Plotting"
    #r.plotTime(N)
    r.labelLeaves(counter=[1])
    with open("test.newick","w") as f:  
        f.write(r.children[0].newick(range(N)))
    # pylab.figure()
    # points=[]
    # r.getDim(0,points)
    # print std(points)
    # pylab.hist(points,N/10)
    # pylab.show()
    

def oldmain():
    N=1000
    C=1.5
    if len(sys.argv)>1:
        N=int(sys.argv[1])
    if len(sys.argv)>2:
        C=float(sys.argv[2])
    if len(sys.argv)>3:
        alpha=float(sys.argv[3])
    if len(sys.argv)>4:
        beta=float(sys.argv[4])
    #test(N,c)
    D=2
    pylab.close('all')
    invA=lambda x: 1.0-exp(-x/C)
    print "Creating ddt"
    r=sampleDDT(N,D,1.0,invA)
    while True:
        pylab.cla()
        r.plot(N)
        #r.diffusePlot(zeros(D),0,1)
        pylab.show()
        raw_input()
        
def oldermain():
    #invA=lambda x: 1.0-exp(-x/C)
    #binarytree=sampleDDT(100,2,1,invA)
    #alpha=1.0
    #nonbinarytree=sampleDDT(100,2,1,invA)
    leaves=[]
    for i in range(4):
        leaves.append(node(0,0))
    left=node(None,None)
    left.children.append(leaves[0])
    left.children.append(leaves[1])
    right=node(None,None)
    right.children.append(leaves[2])
    right.children.append(leaves[3])
    
    binarytree=node(None,None)
    binarytree.children.append(left)
    binarytree.children.append(right)
    
    binarytree.countChildren()
    
    nonbinarytree=node(None,None)
    nonbinarytree.children=leaves
    
    nonbinarytree.countChildren()
    
    alphas=arange(0,2,0.02)
    res=zeros([2,len(alphas)])
    for i in range(len(alphas)):
        alpha=alphas[i]
        #print binarytree.structureProbability() + binarytree.timesProbability()
        res[0,i]=binarytree.structureProbability()
        #print binarytree.timesProbability()
        #print nonbinarytree.structureProbability() + nonbinarytree.timesProbability()
        res[1,i]=nonbinarytree.structureProbability()
        #print nonbinarytree.timesProbability()
        
    pylab.plot(alphas,res[0,:])
    pylab.plot(alphas,res[1,:])
    pylab.xlabel('alpha')
    pylab.ylabel('log probability')
    pylab.show()
    
if __name__=="__main__":
    #oldmain()
    DPMsample(1000)
    exit()
    N=15
    S=10
    C=5/double(S)
    
    #random.seed(1)
    if len(sys.argv)>1:
        N=int(sys.argv[1])
    if len(sys.argv)>2:
        C=float(sys.argv[2])
    if len(sys.argv)>3:
        alpha=float(sys.argv[3])
    if len(sys.argv)>4:
        beta=float(sys.argv[4])
    #test(N,c)
    D=1
    pylab.close('all')
    A=lambda x: C/(1.0+x)
    print "Creating ddt"
    r=sampleNCRP(N,S,A)
    G.write_png('example1_graph.png')