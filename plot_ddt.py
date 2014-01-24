# -*- coding: utf-8 -*-
import sys, os
from numpy import *
import pylab
import logging
pylab.ion()
#basedir="synthetic-fillin/"
#resultsdir="em-cheat/"
basedir="fractal/"
#resultsdir="multiple_trees/"
#basedir="synthetic-final/"
resultsdir="sampling/"
repeats=1
for R in range(1,repeats+1):
    print R
    pylab.cla()
    train=genfromtxt(sys.argv[1],dtype="float")

    #try:
        #density=genfromtxt(basedir+resultsdir+"density%d.txt"%R,delimiter=",",dtype="float")[0:400,0:400]
        #print density.sum()*.01*0.01
    
        ##pylab.contour(arange(-2,2,.01),arange(-2,2,.01),log(density.transpose()))
        #pylab.imshow(log(density.transpose()),cmap=pylab.cm.gray, extent=[-2,2,-2,2], origin='lower')
    #except:
        #print "No density file found"
    #break
    tree=genfromtxt(sys.argv[2],delimiter=" ",dtype="float")
    pylab.plot(train[:,0],train[:,1],'ko',ms=2)
    #break
    #testfn=basedir+"test%d.txt"%R
    #if os.path.exists(testfn):
        #test=genfromtxt(testfn,dtype="float")
        #pylab.scatter(test[:,0],test[:,1],c='g')
    pylab.xlim([-2.0,2.0])
    pylab.ylim([-2.0,2.0])
    toptree=tree[tree[:,0]>=0,:]
    s=toptree.shape[0]
    for i in range(s):
        pylab.plot(toptree[i,1:3],toptree[i,3:5],c="b")

    #toptree2=tree[tree[:,0]>=2,:]
    #pylab.scatter(toptree2[:,1], toptree2[:,3],c='g')

    pylab.figure()
    for i in range(s):
        pylab.plot(toptree[i,5:7],toptree[i,1:3],c="b")
    pylab.draw()
    #raw_input()
pylab.ioff()
    
pylab.show()

    # % n=size(train1,1);
    # % n=1
    # % while (n<size(train4,1))
    # %     scatter(train4(1:n,1),train4(1:n,2),'r')
    # %     input(' ')
    # %     n=n+1;
    # % end
