## NOTE: to run the supplementary material, use Jupyter notebook. See README

import numpy as np
import sys
from scipy.sparse import coo_matrix

def AMFilter(xin, nelx, nely, sensList=[]):
# AMFILTER Applies a virtual additive manufacturing process to a 
#          2D blueprint design input.
#    Possible uses:
#    xi = AMfilter(x,nelx,nely)    design transformation
#    [xi, (df1dx, df2dx,...)] = AMfilter(x, [df1dxi, df2dxi, ...])
#        This includes also the transformation of design sensitivities
#  where
#    x : blueprint design (2D array), 0 <= x(i,j) <= 1
#    xi: printed design (2D array)
#    df1dx, df1dxi etc.:  design sensitivity (2D arrays)

    #INTERNAL SETTINGS
    P = 60.; ep = 1e-4; xi_0 = 0.5; # parameters for smooth max/min functions
    dfx=[]

    nSens = len(sensList); 

    #ORIENTATION
    x = xin.reshape(nelx,nely).T
    xi = np.zeros(x.shape)
    for k in range(0,nSens):
        sensList[k] = sensList[k].reshape(nelx,nely).T
    
    [nely,nelx]=x.shape; 

    #AM FILTER =====================
    Ns = 3
    Q = P + np.log(Ns) / np.log(xi_0); 
    SHIFT = 100 * (np.finfo(float).tiny)**(1/P); # small shift to prevent division by 0
    BACKSHIFT = 0.95 * Ns**(1./Q) * SHIFT**(P/Q)
    Xi=np.zeros(x.shape)
    keep=np.zeros(x.shape) 
    sq=np.zeros(x.shape)
    # baseline: identity
    xi[nely-1,:]=x[nely-1,:].copy(); # copy base row as-is
    cbr = np.zeros([xi.shape[1]+2])
    for i in range((nely-2),-1,-1):
        # compute maxima of current base row
        cbr[1:xi.shape[1]+1] = xi[i+1,:] + SHIFT
        keep[i,:] = np.power(cbr[0:nelx],P) + np.power(cbr[1:(nelx+1)],P) + np.power(cbr[2:],P)
        Xi[i,:] = np.power(keep[i,:],1./Q) - BACKSHIFT
        sq[i,:] = np.sqrt(np.power(x[i,:]-Xi[i,:], 2) + ep)

        # set row above to supported value using smooth minimum:
        xi[i,:] = 0.5*((x[i,:]+Xi[i,:]) - sq[i,:] + np.sqrt(ep))

    #SENSITIVITIES
    if nSens:
        lmbda = np.zeros([nSens, nelx]) 
        for k in range(0,nSens):
            dfx.append(np.zeros(x.shape))

        # from top to base layer:
        for i in range(0,nely-1):
            # smin sensitivity terms
            dsmindx  = .5*(1. - np.divide(x[i,:]-Xi[i,:], sq[i,:]))
            #dsmindXi = .5*(1+(x(i,:)-Xi(i,:))./sq(i,:)); 
            dsmindXi = 1. - dsmindx
            # smax sensitivity terms
            cbr[1:xi.shape[1]+1] = xi[i+1,:] + SHIFT
            dmx = np.zeros([Ns, nelx])
            for j in range(0,Ns):
                dmx[j,:] = (P/Q)*np.multiply(np.power(keep[i,:], 1./Q-1.), np.power(cbr[j:nelx+j], P-1.))
       
            # rearrange data for quick multiplication:
            qj = np.tile(np.array([-1, 0, 1]), (1, nelx))
            qi = np.tile(np.arange(0,nelx),(3,1))
            qi = qi.flatten('F')
            qj = (qj + qi).flatten()
            qs = dmx.flatten('F')
            dsmaxdxi = coo_matrix((qs[2:-1], (qi[2:-1], qj[2:-1])))
            for k in range(0,nSens):
                dfx[k][i,:] = np.multiply(dsmindx, sensList[k][i,:]+lmbda[k,:])
                lmbda[k,:] = np.multiply(sensList[k][i,:]+lmbda[k,:], dsmindXi)*dsmaxdxi

        # base layer:
        i = nely-1
        for k in range(0,nSens):
            dfx[k][i,:] = sensList[k][i,:]+lmbda[k,:]
        
    xi = xi.flatten('F')

    if nSens:
        argsout = []
        for s in range(0,nSens):
            argsout.append(dfx[s].flatten('F'))
        return (xi, argsout)
    else:
        return (xi)

############################################################################
# This Python code is an adaption of the AM filter matlab code             #
# by Matthijs Langelaar. Adaption to Python by Emiel van de Ven.           #
#                                                                          #
# Statement from the orignal code:                                         #
# This Matlab code was written by Matthijs Langelaar,                      #
# Department of Precision and Microsystems Engineering,                    #
# Delft University of Technology, Delft, the Netherlands.                  #
# Please sent your comments to: m.langelaar@tudelft.nl                     #
#                                                                          #
# The code is intended for educational purposes and theoretical details    #
# are discussed in the paper "An additive manufacturing filter for         #
# topology optimization of print-ready designs", M. Langelaar (2016),      #
# Struct Multidisc Optim, DOI: 10.1007/s00158-016-1522-2.                  #
#                                                                          #
# This code is intended for integration in the 88-line topology            #
# optimization code discussed in the paper                                 #
# "Efficient topology optimization in MATLAB using 88 lines of code,       #
# E. Andreassen, A. Clausen, M. Schevenels, B. S. Lazarov and O. Sigmund,  # 
# Struct Multidisc Optim, 2010, Vol 21, pp. 120--127.                      #
#                                                                          #
# Disclaimer:                                                              #
# The author reserves all rights but does not guarantee that the code is   #
# free from errors. Furthermore, the author shall not be liable in any     #
# event caused by the use of the program.                                  #
############################################################################
