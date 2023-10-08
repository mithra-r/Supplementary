# NOTE: to run the supplementary material, use Jupyter notebook. See README

## Disclaimer:                                                              
## The author reserves all rights but does not guarantee that the code is   
## free from errors. Furthermore, the author shall not be liable in any     
## event caused by the use of the program.                                  

# A 200 LINE TOPOLOGY OPTIMIZATION CODE BY NIELS AAGE AND VILLADS EGEDE JOHANSEN, JANUARY 2013
# Updated by Niels Aage February 2016
# Modifications for overhang constraint by Emiel van de Ven 2019
from __future__ import division
import numpy as np

from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve, splu
from matplotlib import colors
import matplotlib.pyplot as plt
from AMFilter import AMFilter
from FPAMFilterWpr import *

def main(opt_controls,am_controls,w_run,w_out):
	nelx=int(opt_controls[0].value)
	nely=int(opt_controls[1].value)
	rmin = float(opt_controls[2].value)
	ft=1
	AMfilterType = am_controls[0].value


	w_volfrac=opt_controls[3]
	w_penal=opt_controls[4]
	w_movelimit=opt_controls[5]
	
	volfrac = w_volfrac.value
	penal = w_penal.value
	cont0 = 1. - float(am_controls[1].value)

	# Max and min stiffness
	Emin=1e-9
	Emax=1.0

	cont = cont0
	if AMfilterType=='front-propagation':
		AMfilt = FPAMFilterCppCreate(nelx, nely, 0, am_controls[2].value, am_controls[3].value)
	elif AMfilterType=='front-propagation-improved':
		AMfilt = FPAMFilterCppCreate(nelx, nely, 1, am_controls[2].value, am_controls[3].value)

	# dofs:
	ndof = 2*(nelx+1)*(nely+1)

	# Allocate design variables (as array), initialize and allocate sens.
	x= volfrac*np.ones(nely*nelx,dtype=float)
	xold=x.copy()
	xPhys=x.copy()

	if AMfilterType=='layer-by-layer':
		xPrint = cont*AMFilter(xPhys,nelx,nely) + (1-cont)*xPhys
	elif AMfilterType=='front-propagation' or AMfilterType=='front-propagation-improved':
		xPrint = cont*FPAMFilterCpp(AMfilt, xPhys) + (1-cont)*xPhys
	elif AMfilterType=='none':
		xPrint = xPhys

	g=0 # must be initialized to use the NGuyen/Paulino OC approach
	dc=np.zeros((nely,nelx), dtype=float)

	# FE: Build the index vectors for the for coo matrix format.
	KE=lk()
	edofMat=np.zeros((nelx*nely,8),dtype=int)
	for elx in range(nelx):
		for ely in range(nely):
			el = ely+elx*nely
			n1=(nely+1)*elx+ely
			n2=(nely+1)*(elx+1)+ely
			edofMat[el,:]=np.array([2*n1+2, 2*n1+3, 2*n2+2, 2*n2+3,2*n2, 2*n2+1, 2*n1, 2*n1+1])
	# Construct the index pointers for the coo format
	iK = np.kron(edofMat,np.ones((8,1))).flatten()
	jK = np.kron(edofMat,np.ones((1,8))).flatten()    

	# Filter: Build (and assemble) the index+data vectors for the coo matrix format
	nfilter=int(nelx*nely*((2*(np.ceil(rmin)-1)+1)**2))
	iH = np.zeros(nfilter)
	jH = np.zeros(nfilter)
	sH = np.zeros(nfilter)
	cc=0
	for i in range(nelx):
		for j in range(nely):
			row=i*nely+j
			kk1=int(np.maximum(i-(np.ceil(rmin)-1),0))
			kk2=int(np.minimum(i+np.ceil(rmin),nelx))
			ll1=int(np.maximum(j-(np.ceil(rmin)-1),0))
			ll2=int(np.minimum(j+np.ceil(rmin),nely))
			for k in range(kk1,kk2):
				for l in range(ll1,ll2):
					col=k*nely+l
					fac=rmin-np.sqrt(((i-k)*(i-k)+(j-l)*(j-l)))
					iH[cc]=row
					jH[cc]=col
					sH[cc]=np.maximum(0.0,fac)
					cc=cc+1
	# Finalize assembly and convert to csc format
	H=coo_matrix((sH,(iH,jH)),shape=(nelx*nely,nelx*nely)).tocsc()	
	Hs=H.sum(1)

	# BC's and support
	dofs=np.arange(2*(nelx+1)*(nely+1))
	fixed=np.union1d(dofs[0:2*(nely+1):2],np.array([2*(nelx+1)*(nely+1)-1]))
	free=np.setdiff1d(dofs,fixed)

	# Solution and RHS vectors
	f=np.zeros((ndof,1))
	u=np.zeros((ndof,1))

	# Set load
	f[1,0]=-1
	# Initialize plot and plot the initial design
	if plt.gcf():
		fig = plt.gcf()
		ax1 = plt.gca()
	else:
		fig,ax1 = plt.subplots()
	im1 = ax1.imshow(-xPrint.reshape((nelx,nely)).T, cmap='gray',\
	interpolation='none',norm=colors.Normalize(vmin=-1,vmax=0))
	fig.tight_layout()
	fig.show()

	loop=0
	change=1
	dv = np.ones(nely*nelx)
	dc = np.ones(nely*nelx)
	ce = np.ones(nely*nelx)
	while change>0.001 and loop<2000 and w_run.value==True:
		loop=loop+1
		# Setup and solve FE problem
		sK=((KE.flatten()[np.newaxis]).T*(Emin+(xPrint)**penal*(Emax-Emin))).flatten(order='F')
		K = coo_matrix((sK,(iK,jK)),shape=(ndof,ndof)).tocsc()
		# Remove constrained dofs from matrix and convert to coo
		K = deleterowcol(K,fixed,fixed)
		#u[free,0] = spsolve(K, f[free,0])
		u[free,0] = splu(K).solve(f[free,0])

		# Objective and sensitivity
		ce[:] = (np.dot(u[edofMat].reshape(nelx*nely,8),KE) * u[edofMat].reshape(nelx*nely,8) ).sum(1)
		obj=( (Emin+xPrint**penal*(Emax-Emin))*ce ).sum()
		dc[:]=(-penal*xPrint**(penal-1)*(Emax-Emin))*ce

		dv[:] = np.ones(nely*nelx)

		
		if AMfilterType=='layer-by-layer':
			[xPrint,(dcp,dvp)] = AMFilter(xPhys,nelx,nely,[dc,dv])
		elif AMfilterType=='front-propagation' or AMfilterType=='front-propagation-improved':
			[dcp,dvp] = FPAMFilterSensCpp(AMfilt, [dc,dv])
		elif AMfilterType=='none':
			dcp = dc
			dvp = dv

		dc = cont*dcp + (1-cont)*dc
		dv = cont*dvp + (1-cont)*dv

		# Sensitivity filtering:
		if ft==0:
			dc[:] = np.asarray((H*(x*dc))[np.newaxis].T/Hs)[:,0] / np.maximum(0.001,x)
		elif ft==1:
			dc[:] = np.asarray(H*(dc[np.newaxis].T/Hs))[:,0]
			dv[:] = np.asarray(H*(dv[np.newaxis].T/Hs))[:,0]

		# Optimality criteria
		xold[:]=x

		dv[dv<0] = 0
		dc[dc>0] = 0

		l1=0.
		l2=1e9
		move=w_movelimit.value
		# reshape to perform vector operations
		xnew=np.zeros(nelx*nely)
		g0 = np.sum(xPrint)
		while (l2-l1)/(l1+l2)>1e-3 and l2>1e-3:
			lmid=0.5*(l2+l1)
			xnew[:]= np.maximum(0.01,np.maximum(x-move,np.minimum(1.0,np.minimum(x+move,x*np.sqrt(-dc/(dv+1e-10)/lmid)))))
			if g0 + np.sum((xnew-x)*dv) > volfrac*nelx*nely :
				l1=lmid
			else:
				l2=lmid

		x[:] = xnew

		if ft==0:   xPhys[:]=xnew
		elif ft==1:	xPhys[:]=np.asarray(H*xnew[np.newaxis].T/Hs)[:,0]
		
		if AMfilterType=='layer-by-layer':
			xPrint = cont*AMFilter(xPhys,nelx,nely) + (1-cont)*xPhys
		elif AMfilterType=='front-propagation' or AMfilterType=='front-propagation-improved':
			xPrint = cont*FPAMFilterCpp(AMfilt, xPhys) + (1-cont)*xPhys
		elif AMfilterType=='none':
			xPrint = xPhys

		if loop>10:
			cont = min(1.0, cont + 0.1)

		# Compute the change by the inf. norm
		change=np.linalg.norm(x.reshape(nelx*nely,1)-xold.reshape(nelx*nely,1),np.inf)

		# Plot to screen
		im1.set_array(-xPrint.reshape((nelx,nely)).T)
		fig.canvas.draw_idle()

		volfrac = w_volfrac.value
		penal = w_penal.value

		# Write iteration history to screen (req. Python 2.6 or newer)
		w_out.value = "it.: {0} , obj.: {1:.3f} Vol.: {2:.3f}, ch.: {3:.3f}".format(\
					loop,obj,(np.sum(xPrint))/(nelx*nely),change)
	# Make sure the plot stays and that the shell remains	
	plt.show()
	if AMfilterType=='front-propagation' or AMfilterType=='front-propagation-improved':
		FPAMFilterCppDelete(AMfilt)
    
#element stiffness matrix
def lk():
	E=1
	nu=0.3
	k=np.array([1/2-nu/6,1/8+nu/8,-1/4-nu/12,-1/8+3*nu/8,-1/4+nu/12,-1/8-nu/8,nu/6,1/8-3*nu/8])
	KE = E/(1-nu**2)*np.array([ [k[0], k[1], k[2], k[3], k[4], k[5], k[6], k[7]],
	[k[1], k[0], k[7], k[6], k[5], k[4], k[3], k[2]],
	[k[2], k[7], k[0], k[5], k[6], k[3], k[4], k[1]],
	[k[3], k[6], k[5], k[0], k[7], k[2], k[1], k[4]],
	[k[4], k[5], k[6], k[7], k[0], k[1], k[2], k[3]],
	[k[5], k[4], k[3], k[2], k[1], k[0], k[7], k[6]],
	[k[6], k[3], k[4], k[1], k[2], k[7], k[0], k[5]],
	[k[7], k[2], k[1], k[4], k[3], k[6], k[5], k[0]] ])
	return (KE)
    
def deleterowcol(A, delrow, delcol):
	# Assumes that matrix is in symmetric csc form !
	m = A.shape[0]
	keep = np.delete (np.arange(0, m), delrow)
	A = A[keep, :]
	keep = np.delete (np.arange(0, m), delcol)
	A = A[:, keep]
	return A    

if __name__== "__main__":
	main(opt_controls, am_controls, w_run, w_out)
