def d_x(x,qq):
    '''
    Return x derivative with fourth order
    Inhomogeneous grid applicable
    
    Parameters
    ----------
    x : numpy.ndarray, float
        x coordinate (1D)
    qq : numpy.ndarray, float
        variable (3D)
    
    Returns
    -------
        qqd : numpy.ndarray, float
            differentiated quantity (3D)
    '''
    import ctypes
    import numpy as np
    import os
    
    mddir = os.path.dirname(__file__)
    libfile = 'fortran_src/derv.so'
    if os.path.isfile(mddir+'/'+libfile):
        lib = np.ctypeslib.load_library('derv.so',mddir+'/fortran_src')
        lib.d_x.argtypes = _derv_array()
        lib.d_x.restype = ctypes.c_void_p
        ix, jx, kx = qq.shape
        
        qqd = np.zeros((ix,jx,kx),dtype=np.float32,order='F')
        lib.d_x(
            np.copy(qq,order='F').astype(np.float32),
            x.astype(np.float32),
            ctypes.byref(ctypes.c_int(ix)),
            ctypes.byref(ctypes.c_int(jx)),
            ctypes.byref(ctypes.c_int(kx)),
            qqd,
        )
        return qqd
    else:
        _warning(mddir,d_x.__name__)

def d_y(y,qq):
    '''
    Return y derivative with fourth order
    Inhomogeneous grid applicable
    
    Parameters
    ----------
    y : numpy.ndarray, float
        y coordinate (1D)
    qq : numpy.ndarray, float
        variable (3D)
    
    Returns
    -------
        qqd : numpy.ndarray, float
            differentiated quantity (3D)
    '''
    import ctypes
    import numpy as np
    import os
    
    mddir = os.path.dirname(__file__)
    libfile = 'fortran_src/derv.so'
    if os.path.isfile(mddir+'/'+libfile):
        lib = np.ctypeslib.load_library('derv.so',mddir+'/fortran_src')
        lib.d_y.argtypes = _derv_array()
        lib.d_y.restype = ctypes.c_void_p
        ix, jx, kx = qq.shape
        
        qqd = np.zeros((ix,jx,kx),dtype=np.float32,order='F')
        lib.d_y(
            np.copy(qq,order='F').astype(np.float32),
            y.astype(np.float32),
            ctypes.byref(ctypes.c_int(ix)),
            ctypes.byref(ctypes.c_int(jx)),
            ctypes.byref(ctypes.c_int(kx)),
            qqd,
        )
        return qqd
    else:
        _warning(mddir,d_y.__name__)

def d_z(z,qq):
    '''
    Return z derivative with fourth order
    Inhomogeneous grid applicable
    
    Parameters
    ----------
    z : numpy.ndarray, float
        z coordinate (1D)
    qq : numpy.ndarray, float
        variable (3D)
    
    Returns
    -------
        qqd : numpy.ndarray, float
            differentiated quantity (3D)
    '''
    import ctypes
    import numpy as np
    import os
    
    mddir = os.path.dirname(__file__)
    libfile = 'fortran_src/derv.so'
    if os.path.isfile(mddir+'/'+libfile):    
        lib = np.ctypeslib.load_library('derv.so',mddir+'/fortran_src')
        lib.d_z.argtypes = _derv_array()
        lib.d_z.restype = ctypes.c_void_p
        ix, jx, kx = qq.shape
        
        qqd = np.zeros((ix,jx,kx),dtype=np.float32,order='F')
        lib.d_z(
            np.copy(qq,order='F').astype(np.float32),
            z.astype(np.float32),
            ctypes.byref(ctypes.c_int(ix)),
            ctypes.byref(ctypes.c_int(jx)),
            ctypes.byref(ctypes.c_int(kx)),
            qqd,
        )
        return qqd
    else:
        _warning(mddir,d_z.__name__)

def  _derv_array():
    """
    Initialize argument types for derv.so
    
    """
    import ctypes
    import numpy as np
    
    argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float32), # qq
        np.ctypeslib.ndpointer(dtype=np.float32), # x
        ctypes.POINTER(ctypes.c_int), # ix
        ctypes.POINTER(ctypes.c_int), # jx
        ctypes.POINTER(ctypes.c_int), # kx
        np.ctypeslib.ndpointer(dtype=np.float32), # qqd
    ]

    return argtypes

def interp(x,y,z,xu,yu,zu,qq):
    '''
    Interpolation of variable qq to xu,yu,zu from x,y,z
    
    Parameters
    ----------
    x : numpy.ndarray, float
        original x coordinate (1D)
    y : numpy.ndarray, float
        original y coordinate (1D)
    z : numpy.ndarray, float
        original z coordinate (1D)
    xu : numpy.ndarray, float
        target x coordinate (1D)
    yu : numpy.ndarray, float
        target y coordinate (1D)
    zu : numpy.ndarray, float
        target z coordinate (1D)
    qq : numpy.ndarray, float
        variable defined on x,y, z(3D)
    
    Returns
    -------
        qqd : numpy.ndarray
            Interpolated quantity (3D)
    '''
    import ctypes
    import numpy as np
    import os
    mddir = os.path.dirname(__file__)
    libfile = 'fortran_src/regrid.so'
    if os.path.isfile(mddir+'/'+libfile):
        lib = np.ctypeslib.load_library('regrid.so',mddir+'/fortran_src')
        lib.interp.argtypes = [
            np.ctypeslib.ndpointer(dtype=np.float64), # x
            np.ctypeslib.ndpointer(dtype=np.float64), # y
            np.ctypeslib.ndpointer(dtype=np.float64), # z
            np.ctypeslib.ndpointer(dtype=np.float64), # xu
            np.ctypeslib.ndpointer(dtype=np.float64), # yu
            np.ctypeslib.ndpointer(dtype=np.float64), # zu
            np.ctypeslib.ndpointer(dtype=np.float64), # qq
            ctypes.POINTER(ctypes.c_int), # ix
            ctypes.POINTER(ctypes.c_int), # jx
            ctypes.POINTER(ctypes.c_int), # kx
            ctypes.POINTER(ctypes.c_int), # ixu
            ctypes.POINTER(ctypes.c_int), # jxu
            ctypes.POINTER(ctypes.c_int), # kxu
            np.ctypeslib.ndpointer(dtype=np.float64), # qqu
        ]
        lib.interp.restype = ctypes.c_void_p
        
        ix,  jx,  kx = len(x) , len(y) , len(z)
        ixu, jxu, kxu= len(xu), len(yu), len(zu)
        qqu = np.empty((ixu,jxu,kxu),dtype=np.float64,order='F')
        lib.interp(
            x,y,z,xu,yu,zu,
            np.copy(qq,order='F').astype(np.float64),
            ctypes.byref(ctypes.c_int(ix)),
            ctypes.byref(ctypes.c_int(jx)),
            ctypes.byref(ctypes.c_int(kx)),
            ctypes.byref(ctypes.c_int(ixu)),
            ctypes.byref(ctypes.c_int(jxu)),
            ctypes.byref(ctypes.c_int(kxu)),
            qqu,
        )
        return qqu
    else:
        _warning(mddir,interp.__name__)
    
def spherical2cartesian(rr,th,ph,qqs,ixc,jxc,kxc):
    '''
    data in spherical geometry is converted to uniform cartesian geometry
    
    Parameters
    ----------
    rr : numpy.ndarray, float
        radius (1D)
    th : numpy.ndarray, float
        colatitude (1D)
    ph : numpy.ndarray, float
        longitude (1D)
    qqs : numpy.ndarray, float
        variable in spherical geometry (3D)
    ixc : int
        No. of grid in converted x coordinate
    jxc : int
        No. of grid in converted y coordinate
    kxc : int
        No. of grid in converted z coordinate
   
    Returns
    -------
    qqc : numpy.ndarray, float
        converted variable (3D, ixc, jxc, kxc)
    xc : numpy.ndarray, float
        converted x coordinate (1D, ixc)
    yc : numpy.ndarray, float
        converted y coordinate (1D, jxc)
    zc : numpy.ndarray, float
        converted z coordinate (1D array, kxc)
        

    Notes
    -----
    The cartesian coordinate is automatically generated with considering
    the computational domain size.

    '''    
    import ctypes
    import numpy as np
    import os    
    mddir = os.path.dirname(__file__)
    libfile = 'fortran_src/regrid.so'
    if os.path.isfile(mddir+'/'+libfile):    
        lib = np.ctypeslib.load_library('geometry_convert.so',mddir+'/fortran_src')    
        lib.spherical2cartesian.argtypes = [
            np.ctypeslib.ndpointer(dtype=np.float64), # qqs
            np.ctypeslib.ndpointer(dtype=np.float64), # rr
            np.ctypeslib.ndpointer(dtype=np.float64), # th
            np.ctypeslib.ndpointer(dtype=np.float64), # ph
            ctypes.POINTER(ctypes.c_int), # ixs
            ctypes.POINTER(ctypes.c_int), # jxs
            ctypes.POINTER(ctypes.c_int), # kxs
            ctypes.POINTER(ctypes.c_int), # ixc
            ctypes.POINTER(ctypes.c_int), # jxc
            ctypes.POINTER(ctypes.c_int), # kxc
            np.ctypeslib.ndpointer(dtype=np.float64), # qqc
            np.ctypeslib.ndpointer(dtype=np.float64), # xc
            np.ctypeslib.ndpointer(dtype=np.float64), # yc
            np.ctypeslib.ndpointer(dtype=np.float64), # zc
        ]
        lib.spherical2cartesian.restype = ctypes.c_void_p
        
        ixs, jxs, kxs = len(rr) , len(th) , len(ph)
        xc = np.empty(ixc,dtype=np.float64)
        yc = np.empty(jxc,dtype=np.float64)
        zc = np.empty(kxc,dtype=np.float64)
        qqc = np.empty((ixc,jxc,kxc),dtype=np.float64,order='F')
        lib.spherical2cartesian(
            np.copy(qqs,order='F').astype(np.float64),
            rr.astype(np.float64),th.astype(np.float64),ph.astype(np.float64),
            ctypes.byref(ctypes.c_int(ixs)),
            ctypes.byref(ctypes.c_int(jxs)),
            ctypes.byref(ctypes.c_int(kxs)),
            ctypes.byref(ctypes.c_int(ixc)),
            ctypes.byref(ctypes.c_int(jxc)),
            ctypes.byref(ctypes.c_int(kxc)),
            qqc,xc,yc,zc,
        )
        return qqc,xc,yc,zc
    else:
        _warning(mddir,spherical2cartesian.__name__)
    
def _warning(mddir,funcname):
    '''
    Warning message for uninstalled function
    
    Parameters
    ----------
    mddir : str
        directory of the module
    
    funcname : str
        function name
        
    Returns
    -------
    None
    '''
    print('This function '+funcname+' has not been installed')
    print('Please make at '+mddir+'/fortran_src')
    return
