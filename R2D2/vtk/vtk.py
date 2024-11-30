######################################################
######################################################
def write_3D(qq,x,y,z,file,name):
    '''
    This function outputs the 3D scalar data in 
    VTK format especially for Paraview

    Parameters:
        qq (float) [ix,jx,kx]: output 3D array
        x (float) [ix]: x coordinate
        y (float) [jx]: y coordinate
        z (float) [kx]: z coordinate
        file (str): File name for output
        name (str): Name of the variable    
    '''
    import os

    ix = qq.shape[0]
    jx = qq.shape[1]
    kx = qq.shape[2]
    
    f = open(file,mode='w')
    f.write('# vtk DataFile Version 3.0\n')
    f.write('vtk_data\n')
    f.write('BINARY\n')
    f.write('DATASET STRUCTURED_POINTS\n')
    f.write('DIMENSIONS '+str(ix)+' '+str(jx)+' '+str(kx)+'\n')
    f.write('ORIGIN '+'{:.8f}'.format(x.min())+' '+'{:.8f}'.format(y.min())+' '+'{:.8f}'.format(z.min())+'\n')
    dx = x[1] - x[0]
    dy = y[1] - y[0]
    dz = z[1] - z[0]
    
    f.write('SPACING '+'{:.8f}'.format(dx)+' '+'{:.8f}'.format(dy)+' '+'{:.8f}'.format(dz)+'\n')
    f.write('POINT_DATA '+str(qq.size)+'\n')
    f.write('SCALARS '+name+' float\n')
    f.write('LOOKUP_TABLE default\n')

    f.close()
    f = open(file,mode='ab')
    f.write(qq.reshape([ix*jx*kx],order='F').astype('>f'))
    #f.write(qqt.astype('>f'))
    f.close()

######################################################
######################################################
def write_optical_surface(qq,height,y,z,file,name):
    '''
    This function outputs the 2D scalar data in 
    VTK format especially for Paraview

    Parameters:
        qq (float) [jx,kx]: output 3D array
        height (float) [jx,kx]: height of optical surface
                                store in self.qt['he']
        y (float) [jx]: y coordinate
        z (float) [kx]: z coordinate
        file (str): File name for output
        name (str): Name of the variable    
    '''
    import os
    import numpy as np

    Y, Z = np.meshgrid(y,z,indexing='ij')
    jx = qq.shape[0]
    kx = qq.shape[1]
    
    xyz = np.zeros((3,jx*kx))
    xyz[0,:] = height.reshape([jx*kx],order='F')
    xyz[1,:] = Y.reshape([jx*kx],order='F')
    xyz[2,:] = Z.reshape([jx*kx],order='F')
    
    
    f = open(file,mode='w')
    f.write('# vtk DataFile Version 3.0\n')
    f.write('vtk_data\n')
    f.write('BINARY\n')
    f.write('DATASET STRUCTURED_GRID\n')
    f.write('DIMENSIONS 1 '+str(jx)+' '+str(kx)+'\n')
    f.write('POINTS '+str(qq.size)+' float\n')
    f.close()

    f = open(file,mode='ab')
    f.write(xyz.reshape([3*jx*kx],order='F').astype('>f'))
    f.close()
    
    f = open(file,mode='a')
    
    f.write('POINT_DATA '+str(qq.size)+'\n')
    f.write('SCALARS '+name+' float\n')
    f.write('LOOKUP_TABLE default\n')

    f.close()
    
    f = open(file,mode='ab')
    f.write(qq.reshape([jx*kx],order='F').astype('>f'))
    f.close()
    
######################################################
######################################################
def write_3D_vector(qx,qy,qz,x,y,z,file,name):
    '''
    This function outputs the 3D vector data in 
    VTK format especially for Paraview

    Parameters:
        qx (float) [ix,jx,kx]: x-component vector
        qy (float) [ix,jx,kx]: y-component vector
        qz (float) [ix,jx,kx]: z-component vector
        x (float) [ix]: x coordinate
        y (float) [jx]: y coordinate
        z (float) [kx]: z coordinate
        file (str): File name for output
        name (str): Name of the variable    
    '''
    import os
    import numpy as np

    ix = qx.shape[0]
    jx = qx.shape[1]
    kx = qx.shape[2]
    
    f = open(file,mode='w')
    f.write('# vtk DataFile Version 3.0\n')
    f.write('vtk_data\n')
    f.write('BINARY\n')
    f.write('DATASET STRUCTURED_POINTS\n')
    f.write('DIMENSIONS '+str(ix)+' '+str(jx)+' '+str(kx)+'\n')
    f.write('ORIGIN '+'{:.8f}'.format(x.min())+' '+'{:.8f}'.format(y.min())+' '+'{:.8f}'.format(z.min())+'\n')
    dx = x[1] - x[0]
    dy = y[1] - y[0]
    dz = z[1] - z[0]
    
    f.write('SPACING '+'{:.8f}'.format(dx)+' '+'{:.8f}'.format(dy)+' '+'{:.8f}'.format(dz)+'\n')
    f.write('POINT_DATA '+str(qx.size)+'\n')
    f.write('VECTORS '+name+' float\n')
    #f.write('SCALARS '+name+' float\n')
    #f.write('LOOKUP_TABLE default\n')
    
    f.close()
    
    f = open(file,mode='ab')
    f.write(np.stack([qx,qy,qz],0).reshape([3*ix*jx*kx],order='F').astype('>f'))        
    #f.write(qx.reshape([ix*jx*kx],order='F').astype('>f'))        
    f.close()

    
