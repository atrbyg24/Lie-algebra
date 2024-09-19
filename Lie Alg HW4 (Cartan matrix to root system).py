#!/usr/bin/env python
# coding: utf-8

# In[2]:



from sage.all import *

#Computes the height of a root
def height(rt):
    ht = 0;
    n = rt.ncols();
    for x in range(n):
        ht = ht + rt[0][x];
    return ht


"""Compute the inner product <beta,alpha_j>
beta is an arbitrary positive root
alpha_j is a simple root"""
def prod(rt,j,Phi,M,B):
    s = 0;
    n = rt.ncols();
    for i in range(n):
        s = s + rt[0][i]*(M[i][j])
    return s;


"""Compute the start of a root string
Find r such that the alpha_j root string through beta is 
beta - r*alpha_j,...,beta + q*alpha_j"""
def root_r(rt,j,Phi,M,B):
    r = 0;
    while(rt in Phi):
        rt = rt - B[j];
        r = r+1;
    return r-1


"""Compute the start of a root string
Find q such that the alpha_j root string through beta is 
beta - r*alpha_j,...,beta + q*alpha_j """    
def root_q(rt,j,Phi,M,B):
    return root_r(rt,j,Phi,M,B)-prod(rt,j,Phi,M,B)

"""Function that takes a Cartan matrix M as the input
Returns the set of positive roots in the root system
Uses root string algorithm as described in 11.1 of Humphreys"""
def CartanToRoot(M):
    n = M.nrows();
    
    #Create space of 1 x n matrices
    MS = MatrixSpace(ZZ,1,n)
    
    """"Get standard basis for MS
    Basis elements of the form [0,0,...0,1,0,...0]
    Basis elements correspond to the simple roots in the root system"""
    B = MS.basis()

    #Create a list for all the positive roots
    Phi = [];
    for v in B:
        Phi.append(v)

    """Compute the alpha_j root string through alpha_i
    Start with height 1 roots, i.e. simple roots
    Given an alpha_j root string through alpha_i
    alpha_i - r*alpha_j,...,alpha_i + q*alpha_j
    r - q = <alpha_i,alpha_j> """
    maxHeight = 1;
    for i in range(n):
        for j in range(n):
            if (M[i][j]<0):
                q = -1*M[i][j];
                if (q + 1 > maxHeight):
                    maxHeight = q + 1;
                for k in range(q):
                    newroot = Phi[i]+Phi[j];
                    if (newroot not in Phi):
                        Phi.append(newroot);
                    
    """Iterate the algorithm starting with height 2 roots
    Algorithm terminates when no new roots can be found
    i.e. all possible root strings have been checked
    maxHeight and curHeight variables determine when algorithm terminates"""
    curHeight = 2;
    while (curHeight <= maxHeight):
        for rt in Phi:
            if (height(rt)==curHeight):
                for j in range(n):
                    r = root_r(rt,j,Phi,M,B);
                    s = prod(rt,j,Phi,M,B)
                    q = root_q(rt,j,Phi,M,B);
                    if (q+height(rt) > maxHeight):
                        maxHeight = q + height(rt);
                    for k in range(q):
                        newroot = rt + B[j];
                    if (newroot not in Phi):
                        Phi.append(newroot);
        curHeight = curHeight + 1;
        
    #Return the positive roots in the root system
    print("The positive roots are ")
    print(Phi)
    print("The number of positive roots is")
    print(len(Phi))


# In[2]:


#Root system for C_3
CartanToRoot(Matrix([[2,-1,0],[-1,2,-1],[0,-2,2]]))


# In[3]:


#Root system for A_4
CartanToRoot(Matrix([[2,-1,0,0],[-1,2,-1,0],[0,-1,2,-1],[0,0,-1,2]]))


# In[4]:


#Root system for G_2
CartanToRoot(Matrix([[2,-1],[-3,2]]))


# In[5]:


#Root system for F_4
CartanToRoot(Matrix([[2,-1,0,0],[-1,2,-2,0],[0,-1,2,-1],[0,0,-1,2]]))


# In[6]:


#Root system for E_6
CartanToRoot(Matrix([[2,0,-1,0,0,0],[0,2,0,-1,0,0],[-1,0,2,-1,0,0],[0,-1,-1,2,-1,0],[0,0,0,-1,2,-1],[0,0,0,0,-1,2]]))


# In[231]:


#Root system for D_4
CartanToRoot(Matrix([[2,-1,0,0],[-1,2,-1,1],[0,-1,2,0],[0,-1,0,2]]))


# In[232]:


#Root system for B_2
CartanToRoot(Matrix([[2,-2],[-1,2]]))


# In[ ]:




