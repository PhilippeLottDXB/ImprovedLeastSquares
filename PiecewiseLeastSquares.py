#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 22 21:04:51 2020

@author: root
"""

#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Sep  8 17:26:08 2018

@author: philippe
"""

import numpy as np
import numpy.polynomial as poly
from copy import deepcopy as dc

import matplotlib.pylab as plt

def lsfit(x,y,order):
    xx  =[]
    #print(x)
    for i in range(order+1):
        xx.append(x**i)
    #print(len(xx))
    #print(len(y),len(x))
    A = np.vstack(xx).T
    m = np.linalg.lstsq(A, y)[0]
    #print(m)
    #print(c)


class derivative_based_least_squares:
    def __init__(self, x, y, dxdt, dydt, order):
        return None





class ConstrainedLeastSquares:
    def __init__(self,x,y,Constraints,order,Continuous=True):
        self.x = x
        self.y = y
        self.order = order
        self.Constraints = Constraints
        ## Assumption : Constraints are all within the system, you cannot put
        ## a constraint outside
        self.numberOfCoefficients = (1+len(Constraints))*(order+1)
        self.powers = np.arange(0,order+1,1)
        self.pCoeffs = np.ones(order+1)
        self.pCoeffs2 = (order+1)*[0]
        #coefficients = np.ones(self.numberOfCoefficients)
        ConstrainedLeastSquares.fillMatrix(self)
        ConstrainedLeastSquares.resolve_flat(self)
        ConstrainedLeastSquares.createLSSystem(self)
        ConstrainedLeastSquares.replaceInMat(self)
        if Continuous:
            polys = []
            for r in range(1+len(Constraints)):
                #print(self.sol4[r*(self.order+1):(1+r)*(self.order+1)])
                newPoly = poly.polynomial.Polynomial(self.sol4[r*(self.order+1):(1+r)*(self.order+1)])
                newPoly = poly.polynomial.Polynomial(self.sol5[r*(self.order+1):(1+r)*(self.order+1)])

                #print(newPoly.coef)
                polys.append(newPoly)
        else:
            polys = []
            for r in range(1+len(Constraints)):
                newPoly = poly.polynomial.Polynomial(self.sol[r*(self.order+1):(1+r)*(self.order+1)])
                polys.append(newPoly)   
        
                
        self.polys = polys
        
    def getIndex(self,x):
        cc= 0
        for c in self.Constraints:
            if x<c[0]:
                return cc
            cc+=1
        return cc
    
    def __call__(self,x):
        cc = ConstrainedLeastSquares.getIndex(self,x)
        return self.polys[cc](x)
        
    def replaceInMat(self):
        s = sum(self.keeper)
        lst0 = []
        lst1 = []
        for j,e in enumerate(self.extractor):
            if e == False:
                lst0.append(j)
            else:
                lst1.append(j)                

        self.Mat3 = np.zeros((s,s))
        ct= 0
        for i,k in enumerate(self.keeper):
            if k:
                self.Mat3[ct]=self.Mat2[i][self.keeper]
                ct+=1
        
        ct = 0

        self.T3 = self.T2[self.keeper] 
        print(len(self.T3))
        print(len(self.cstes))
        for i,m in enumerate(self.Mat2):
            if self.keeper[i]:
                for j,li in enumerate(lst1):
                    self.Mat3[ct] = self.Mat3[ct]+self.Mat2[i,li]*self.m1[j]#-self.cstes[j]
                    self.T3[ct]=self.T3[ct]+self.Mat2[i,li]*self.cstes[j]
                #print(ct)
                #self.T3[ct]=self.T3[ct]+self.cstes[ct]
                ct+=1
        #print self.Mat3
        self.sol4 = np.zeros_like(self.T2)
        #print self.Mat3.shape, self.T3.shape
        self.sol3 = np.linalg.solve(self.Mat3,self.T3)
        self.sol4[self.keeper]=self.sol3
        temp = np.ones_like(self.m2[0])
        for i,k in enumerate(self.m1):
            temp[i]=sum(k*self.sol3)+self.cstes[i]
        self.sol4[self.extractor]=temp
        
        
    def createLSSystem(self):
        xs = []
        ys = []
        b_ = self.x<= self.Constraints[0][0]
        xs.append(np.asarray(self.x[b_]))
        ys.append(np.asarray(self.y[b_]))
        for c1,c2 in zip(self.Constraints[:-1],self.Constraints[1:]):
            b1_ = self.x>= c1[0]
            b2_ = self.x<= c2[0]
            b_  = np.logical_and(b1_,b2_)
            xs.append(np.asarray(self.x[b_]))
            ys.append(np.asarray(self.y[b_]))       
        b_ = self.x>= self.Constraints[-1][0]
        xs.append(np.asarray(self.x[b_]))
        ys.append(np.asarray(self.y[b_]))
        self.xs = xs
        self.ys = ys
        #print(len(xs))
        #print(len(ys))
        TS = len(self.xs)*(1+self.order)
        #print(TS)
        #print(self.NConstraints)
        Mat2 = np.zeros((TS,TS))
        T2   = np.zeros(TS)
        print("TS",TS,"CSTEST",len(self.cstes))
        for i in range(len(self.xs)):
            #lsfit(self.xs[i],self.ys[i],self.order)
            for j in range(self.order+1):
                #print(xs[i])
                #print(ys[i])
                a1 = sum((xs[i]**j)*ys[i])
                for k in range(self.order+1):
                    #selfprint(i,j,k)
                    a2 = sum((xs[i]**j)*(xs[i]**k))
                    Mat2[i*(self.order+1)+j][i*(self.order+1)+k]=a2
                T2[i*(self.order+1)+j] = a1
                
        self.sol = np.linalg.solve(Mat2,T2)
        self.Mat2 = Mat2
        self.T2  = T2
        self.Mat5 = dc(Mat2)
        self.Mat5[self.extractor,:]=self.Matrix
        self.T5 = dc(T2)
        self.T5[self.extractor]=self.cstes
        self.sol5 = np.linalg.solve(self.Mat5,self.T5)
     
    def fillMatrix(self):
        """
        This tool takes the constraints and generates the solving materials for this 
        """
        Matrix = []
        cstes  = []
        extractor = [True]+(self.order)*[False]#+[True]
        zeros = (self.order+1)*[0.0]
        nc = 0
        for i,c in enumerate(self.Constraints):
            print(i)
            lc1 = len(c[1:])
            #extractor+=(self.order+1-lc1)*[False]+lc1*[True]
            xx = c[0]
            lc1 = 0
            lc2 = 0 
            for j,c2 in enumerate(c[1:]):
                lc2+=1
                pre = i * self.pCoeffs2
                post = (len(self.Constraints)-1-i)*self.pCoeffs2
                pows = dc(self.powers)
                coes = dc(self.pCoeffs)
                for k in range(j):
                    coes = coes*pows
                    pows -=1
                if  type(c2) is bool:#True:
                    if c2:
                        nc+=1
                        lc1+=1
                        cstes.append(0.0)
                        rhs = (coes*xx**pows).tolist()
                        lhs = (-coes*xx**pows).tolist()
                        line = pre+rhs+lhs+post
                        Matrix.append(line)

                else:# c2 is not True:
                    print(c2)
                    cstes.append(c2)
                    nc+=1
                    lc1+=1
                    rhs = (coes*xx**pows).tolist()
                    lhs = zeros
                    line = pre+rhs+lhs+post
                    Matrix.append(line)
                    
                    nc+=1
                    lc1+=1                    
                    cstes.append(c2)
                    rhs = zeros
                    lhs = (coes*xx**pows).tolist()
                    line = pre+rhs+lhs+post
                    Matrix.append(line)
                    
                    
            #print(len(Matrix[-1]))
            extractor+=(0)*[True]+(self.order+1-lc1)*[False]+lc1*[True]
        
        #print("Num Constraints",nc)
        extractor[-1]=False
        self.NConstraints = nc
        #print("tests")
        #for e in Matrix:
        #    print(len(e))
        self.Matrix=np.asarray(Matrix)
        #print("Shape Maztix",self.Matrix.shape)
        self.extractor=np.asarray(extractor)
        self.keeper   =np.logical_not(self.extractor)
        self.cstes = np.asarray(cstes)
        
        #print("Shape Keeper",self.keeper.shape)
        #print("Shape Matrix",self.Matrix.shape)
        self.m1 = self.Matrix[:,self.keeper]
        self.m2 = -1*self.Matrix[:,self.extractor]
        #self.m_test = self.Matrix[self.extractor,:]

    def resolve_flat(self):
        nc = self.NConstraints
        for j in np.arange(0,nc,1):
            for i in np.arange(0,nc,1):
                if i!=j:
                    if self.m2[i,j]!=0:
                        koeff = self.m2[j,j]/self.m2[i,j]
                        koeff=1/koeff
                        self.m1[i,:] =self.m1[i,:]-koeff*self.m1[j,:]
                        self.m2[i,:] =self.m2[i,:]-koeff*self.m2[j,:]
                        self.cstes[i]=self.cstes[i]-koeff*self.cstes[j]
            self.m1[j,:]=self.m1[j,:]/self.m2[j,j]
            self.cstes[j]=self.cstes[j]/self.m2[j,j]            
            self.m2[j,:]=self.m2[j,:]/self.m2[j,j]
        #print(self.cstes)
            

nn = 201
x = np.linspace(0,8,nn)
y = x*x-x+0.5*x*x*x*x*x+np.sin(x)+np.random.random_sample(nn)
y = np.sin(0.5*np.pi*x)#+np.random.random_sample(nn)*1.5
y = np.sqrt(16.1-(x-4)**2)*.5+np.random.random_sample(nn)*1.9+y-2.0
settings = [[1.0,0,1,2,3],[2.0,0,1,2,3],[3.0,0,1,2,3],[4.0,0,1,2,3],[5.0,0,1,2,3],[6.0,0,1,2,3],[7.0,0,1,2,3]]
settings = [[1.5,0,1,2],[2.0,0,1,2],[3.5,0,1,2],[4.0,0,1,2],[5.0,0,1,2],[6.0,0,1,2],[7.0,0,1,2]]
settings = [[1.5,True,True,True],[2.0,True,True,True],[3.5,True,True,True],[4.0,True,True,True],[5.0,True,True,True],[6.0,True,True,True],[7.0,True,True,True]]
#settings = [[1.5,True,True,True],[2.0,0.25,True,True],[3.5,-0.8,True,True],[4.0,True,True,True],[5.0,True,True,True],[6.0,True,True,True],[7.0,True,True,True]]
settings = [[1.5,True,True],[2.0,True,True],[3.5,True,True],
            [4.0,True,True],[5.0,True,True],[6.0,True,True],
            [7.0,False,True]]
#r=ConstrainedLeastSquares(x,y,settings,4)
print("TEST 1 successful!")
if False:
    settings = [[1.0,0.0,True],[2.0,1.0,True],[3.,0.0,True],
                [4.0,1.0,True],[5.0,0.0,True],[6.0,1.0,True],
                [7.0,-1.0,True]]
elif False:
    settings = [[1.0,1.0],[2.0,0.0],[3.,-1],
                [4.0,0],[5.0,1],[6.0,0.0],
                [7.0,-1.0]]
elif False:
    settings = [[1.0,1.0,True],[2.0,0.0,True],[3.,-1,True],
                [4.0,0,True],[5.0,1,True],[6.0,0.0,True],
                [7.0,-1.0,True]]
elif False:
    settings = [[1.0,True,True],[2.0,True,True],[3.,0.0,False],
                [4.0,True,True],[5.0,True,True],[6.0,True,True],
                [7.0,1.5,True]] 
elif False:
    settings = [[1.0,True,True,True],[2.0,True,True,True],[3.,True,True,True],
                [4.0,True,True,True],[5.0,True,True,True],[6.0,True,True,True],
                [7.0,True,True,True]]    
elif True:
    settings = [[1.0,True,True,True],[2.0,True,True,True],[3.,True,True,True],
                [4.0,True,True,True],[5.0,True,True,True],[6.0,True,True,True],
                [7.0,True,True,True]]    
    
elif False:
    settings = [[1.0,False,False,True],[2.0,0.0,True,True],[3.,True,True,True],
                [4.0,False,False,True],[5.0,False,False,True],[6.0,False,False,True],
                [7.0,False,False,True]]     
else :# False:
    settings = [[1.0,True,True,True,True],[2.0,True,True,True,True],[3.,True,True,True,True],
                [4.0,True,True,True,True],[5.0,True,True,True,True],[6.0,True,True,True,True],
                [7.0,True,True,True,True]]  

#settings = [[1.0,True,True,True],[4.0,True,True,True],[7.0,True,True,True]] 
#settings = [[.1,1.5],[4.0,-2.2],[7.0,1.5]] 
#settings = [[1.0,True],[4.0,True],[7.0,True]] 

r=ConstrainedLeastSquares(x,y,settings,3)
#settings = [[1.0,0,1],[2.0,0,1],[3.0,0,1],[4.0,0,1],[5.0,0,1],[6.0,0,1],[7.0,0,1]]
#settings = [[1.0,0],[2.0,0],[3.0,0],[4.0,0],[5.0,0],[6.0,0],[7.0,0]]
ranges = []
ct = 0
start = True
for s1,s2 in zip([[0.0]]+settings,settings+[[8.0]]):
    #print(s1,s2)
    ranges.append(np.linspace(s1[0],s2[0],40)[:])
    if start:
        print(ct,r(s1[0]))
        ct+=1
        print(ct,r(s2[0]))
        ct+=1
        start=False
    else:
        print(ct,r(s2[0]))
        ct+=1

#r=ConstrainedLeastSquares(x,y,
#                          [[0.5,0],[1.0,0],[1.5,0],[2.5,0],[3.0,0],[3.5,0]],3)

#r=ConstrainedLeastSquares(x,y,settings,3)
plt.figure(figsize=(11,6))
plt.plot(x,y,'k+',label="raw data")
xs = []
for i,vals in enumerate(ranges):
    ee = []
    for j in vals:
        ee.append(r(j))
    plt.plot(vals,ee,"-",label="result %d"%(1+i))
plt.ylim([-3.5,3.5])
plt.legend(loc="best")
plt.show()


