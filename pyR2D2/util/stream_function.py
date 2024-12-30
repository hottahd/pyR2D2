def stream_function(vrr,vth,rr,th):
    '''
    Calculate stream function from vector field in 2D pherical geometry.
    The vector field must satisfy div v = 0 condition.
    
    Parameters
    ----------
    vrr: numpy.ndarray, float
        Radial velocity field
    vth: numpy.ndarray, float
        Colatitudinal velocity field
    rr: numpy.ndarray, float
        Radial grid, float
    th: numpy.ndarray
        Colatitudinal grid
        
    Returns
    -------
    ph: numpy.ndarray, float
        Stream function
        
    Notes
    -----
    rr and th must be uniform grids.
    '''
    import numpy as np
    drr = rr[1] - rr[0]
    dth = th[1] - th[0]

    ix = vrr.shape[0]
    jx = vrr.shape[1]

    ixg = ix + 2
    jxg = jx + 2

    # grid with margin
    rrg = np.zeros((ixg))
    thg = np.zeros((jxg))
    
    # vector field with margin
    vrrg = np.zeros((ixg, jxg))
    vthg = np.zeros((ixg, jxg))

    rrg[1:ixg-1] = rr
    thg[1:jxg-1] = th

    vrrg[1:ixg-1,1:jxg-1] = vrr
    vthg[1:ixg-1,1:jxg-1] = vth

    rrg[0] = rrg[1] - (rrg[2]-rrg[1])
    rrg[ixg-1] = rrg[ixg-2] + (rrg[ixg-2]-rrg[ixg-3])

    thg[0] = thg[1] - (thg[2]-thg[1])
    thg[jxg-1] = thg[jxg-2] + (thg[jxg-2]-thg[jxg-3])

    RRg, THg = np.meshgrid(rrg,thg,indexing='ij')

    # boundary condition (will be modified)
    ## Radial boundary condition
    vrrg[0,:] = -vrrg[1,:]
    vrrg[ixg-1,:] = -vrrg[ixg-2,:]

    vthg[0,:] = vthg[1,:]
    vthg[ixg-1,:] = vthg[ixg-2,:]

    ## Colatitudinal boundary condition
    vrrg[:,0] = vrrg[:,1]
    vrrg[:,jxg-1] = vrrg[:,jxg-2]

    vthg[:,0] = -vthg[:,1]
    vthg[:,jxg-1] = -vthg[:,jxg-2]

    # vorticity
    ffg = np.zeros((ixg, jxg))

    rvthg = RRg*vthg
    ## calculating vorticity at i,j
    # ffg at the boundary is zero thus we do not set the value additionally
    for j in range(1,jxg-1):
        for i in range(1,ixg-1):
            ffg[i,j] = +0.5*(rvthg[i+1,j] +rvthg[i-1,j])/drr/rrg[i] \
                       -0.5*( vrrg[i,j+1] + vrrg[i,j-1])/dth/rrg[i]
                        
    phm = np.zeros((ixg,jxg))
    phn = np.zeros((ixg,jxg))

    ## initial condition
    # initial guess is based on simple integration
    # calculation based on vth
    for j in range(0,jxg-1):
        for i in range(1,ixg-1):
            phn[i,j] = (rrg[i-1]*phn[i-1,j] - 0.5*(rvthg[i,j] + rvthg[i-1,j])*drr)/rrg[i]
        
    phmx = np.zeros((ixg,jxg))
    phmy = np.zeros((ixg,jxg))
    for m in range(0,5):
        phm = phn.copy()

        # boundary condition (will be modified)
        phm[0,:] = 0.0
        phm[ixg-1,:] = 0.0
        
        phm[:,0] = 0.0
        phm[:,jxg-1] = 0.0
    
        phn[0,:] = 0.0
        phn[ixg-1,:] = 0.0
        
        phn[:,0] = 0.0
        phn[:,jxg-1] = 0.0
        
        dt = drr**2/4

        for j in range(1,jxg):
            for i in range(1,ixg):
                rrm = 0.5*(rrg[i] + rrg[i-1])
                thm = 0.5*(thg[j] + thg[j-1])
                phmx[i,j] = rrm**2     *(phm[i,j] - phm[i-1,j])/drr
                phmy[i,j] = np.sin(thm)*(phm[i,j] - phm[i,j-1])/dth
            
        for j in range(1,jxg-1):
            for i in range(1,ixg-1):
                phmxx = (phmx[i+1,j] - phmx[i,j])/drr/rrg[i]**2
                phmyy = (phmy[i,j+1] - phmy[i,j])/dth/rrg[i]**2/np.sin(thg[j])
                laplace_phm = phmxx + phmyy

                phn[i,j] = phm[i,j] \
                    + (laplace_phm - phm[i,j]/rrg[i]**2/np.sin(thg[j]) + ffg[i,j])*dt
            
    return phn[1:ixg-1,1:jxg-1]            
