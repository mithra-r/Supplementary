## NOTE: to run the supplementary material, use Jupyter notebook. See README

## Wrapper for FPAMFilter.cpp such that its functions can be called from python
##
## Part of supplementary material of "Overhang control in topology optimization: 
## a comparison of continuous front propagation-based and discrete layer-by-layer
## overhang control", E. van de Ven, R. Maas, C. Ayas, M. Langelaar, F. van Keulen,
## 2019
##
## Disclaimer:                                                              
## The author reserves all rights but does not guarantee that the code is   
## free from errors. Furthermore, the author shall not be liable in any     
## event caused by the use of the program.                                  
##
## Code by Emiel van de Ven, 2020
## emiel@emielvandeven.nl

from ctypes import c_int,c_double,c_void_p,cdll,POINTER
import numpy as np
import platform

if "Windows" in platform.system(): 
    from ctypes import WinDLL
    try:
        lib = WinDLL ("FPAMFilter64.dll")
    except:
        print("Could not load 64bit dll, trying 32 bit...")
        try:
            lib = WinDLL ("FPAMFilter32.dll")
        except:
            print("ERROR: Could not load 32bit dll. Try compiling FPAMFilter.cpp yourself.")
elif "Darwin" in platform.system():
    try:
        lib = cdll.LoadLibrary('./libFPAMFilterOSX.so')
    except:
        print("ERROR: Could not load libFPAMFilterOSX.so. Try compiling FPAMFilter.cpp yourself.")
else:
    try:
        lib = cdll.LoadLibrary('./libFPAMFilterUNIX.so')
    except:
        print("ERROR: Could not load libFPAMFilterUNIX.so. Try compiling FPAMFilter.cpp yourself.")

# define input and output types for functions in libfpam.so
lib.FPAMFilter_new.argtypes = c_int,c_int,c_int,c_double,c_double
lib.FPAMFilter_new.restype = c_void_p
lib.FPAMFilter_delete.argtype = c_void_p
lib.FPAMFilter_delete.restype = None
lib.evaluate.argtypes = c_void_p,POINTER(c_double),POINTER(c_double)
lib.evaluate.restype = None
lib.sens.argtypes = c_void_p,POINTER(c_double),POINTER(c_double)
lib.sens.restype = None

# create FPAMFiler object
def FPAMFilterCppCreate(nelx, nely, speedFunctionType, kappa, v_void):
    fpamFilter = c_void_p(lib.FPAMFilter_new(nelx,nely, speedFunctionType,kappa,v_void))
    return fpamFilter

# delete FPAMFiler object
def FPAMFilterCppDelete(fpamFilter):
    lib.FPAMFilter_delete(fpamFilter)

# evaluate AMFilter
def FPAMFilterCpp(fpamFilter, xinPy):
    # cast numpy arrays to correct type
    xin = (c_double * len(xinPy))(*xinPy)
    xout = (c_double * len(xinPy))()

    lib.evaluate(fpamFilter,xin, xout) 

    oh_filtered = np.array(xout)
    return oh_filtered

# evaluate sensitivities
def FPAMFilterSensCpp(fpamFilter, sensList):
    sensOut=[]
    # loop over sensitivities
    for k in range(0,len(sensList)):
        # cast numpy arrays to correct type
        xin = (c_double * len(sensList[k]))(*sensList[k])
        xout = (c_double * len(sensList[k]))()

        lib.sens(fpamFilter, xin, xout) 

        sensOut.append(np.array(xout))

    return sensOut
