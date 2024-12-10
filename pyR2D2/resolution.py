######################################################
######################################################
def gen_coord(xmax,xmin,ix,margin):
    '''
    This fucntion defines uniform geometry

    Parameters:
        xmax (float): location of upper boundary
        xmin (float): location of lower boundary
        ix (int): number of grid without margin
        margin (int): number of margin

    Return:
        x (float) [ix + 2*margin]: generated geometry

    '''
    import numpy as np
    dx = (xmax - xmin)/ix
    x = np.zeros(ix+margin*2)
    x[0] = xmin - (margin - 0.5)*dx
    for i in range(1,ix+2*margin):
        x[i] = x[i-1] + dx
    #x = np.arange(xmin - (margin - 0.5)*dx,xmax + (margin + 0.5)*dx,dx)

    return x

######################################################
######################################################
def gen_coord_ununiform_top(xmax,xmin,ix,margin,dx00,ix_ununi):
    '''
    This function defines ununiform geometry for vertical geometry
    The fine grid is concentrated around the top boundary

    Parameters: 
        xmax (float): location of upper boundary
        xmin (float): location of lower boundary
        ix (int): number of grid without margin
        margin (int): number of margin
        dx00 (float): grid spacing in uniform grid region
        ix_uniuni (int): number of uniform grid

    Return:
        x (float) [ix + 2*margin]: generated geometry
   
    '''
    import numpy as np
    
    xrange = xmax - xmin
    xrange0 = dx00*ix_ununi
    xrange1 = xrange - xrange0
    nxx = ix - ix_ununi

    ixg = ix + 2*margin

    fdx = 2*(xrange1 - dx00*nxx)/(nxx-4)/nxx
    x = np.zeros(ixg)
    x[ixg - margin - 1] = xmax - 0.5*dx00
    for i in range(ixg - margin,ixg):
        x[i] = x[i-1] + dx00

    for i in range(ixg - margin - 2, ixg - margin - ix_ununi - 3, -1):
        x[i] = x[i+1] - dx00

    dx11 = dx00
    for i in range(ixg - margin - ix_ununi-3,3,-1):
        x[i] = x[i+1] - dx11
        dx11 = dx11 + fdx

    for i in range(3,-1,-1):
        x[i] = x[i+1] - dx11
    
    return x

######################################################
######################################################
def gen_coord_ununiform_flex(xmax,xmin,ix,margin,dx_fine,ix_fine,xc_fine):
    '''
    This function defines ununiform geometry
    The fine grid is concentrated around x = xc

    Parameters:
        xmax (float): location of upper boundary
        xmin (float): location of lower boundary
        ix (int): number of grid without margin
        margin (int): number of margin
        dx_fine (float): grid spacing in fine grid
        ix_fine (int): number of fine grid
        xc_fine (float): center position of fine grid

    Return:
        x (float) [ix + 2*margin]: generated geometry        
    '''
    import numpy as np
    
    ixg = ix + 2*margin
    x = np.zeros(ixg)

    # number of low resolution grid
    ix_low = ix - ix_fine
    
    # grid spacing in low resolution grid
    dx_low = (xmax - xmin - dx_fine*ix_fine)/ix_low
    
    if xc_fine - xmin <= dx_fine*ix_fine/2 : # if fine grid cross the boundary
        if((xc_fine - xmin)//dx_fine*dx_fine == xc_fine - xmin):
            xc_fine = (((xc_fine - xmin)//dx_fine)  )*dx_fine + xmin
        else:
            xc_fine = (((xc_fine - xmin)//dx_fine)+1)*dx_fine + xmin
        # number of shift
        ix_shift = np.int((xc_fine - xmin)//dx_fine)
    else: # if coarse grid cross the boundary
        xct = xc_fine - dx_fine*ix_fine/2
        if (xct - xmin)//dx_low*dx_low == xct - xmin:
            xct = ((xct - xmin)//dx_low  )*dx_low + xmin
        else:
            xct = ((xct - xmin)//dx_low+1)*dx_low + xmin
        xc_fine = xct + dx_fine*ix_fine/2
        # number of shift
        ix_shift = np.int((xct - xmin)//dx_low) + ix_fine//2

    # Grid construction start from center of fine grid
    x[margin + ix_shift] = xc_fine + 0.5*dx_fine
    for i in range(1,ix):
        ii = i + ix_shift
        # if the grid exceeds the top boundary
        # return to the bottom boundary
        if i + ix_shift > ix - 1:
            ii = i + ix_shift - ix

        # return procedure from top to bottom boundary
        if ii == 0:
            x[margin - 1] = x[margin + ix - 1] - (xmax - xmin)
            
        if i == ix_fine/2 or i == ix_fine/2 + ix_low: # for fine and coarse grid boundary
            x[margin + ii] = x[margin + ii - 1] + 0.5*(dx_fine + dx_low)
        elif i <= ix_fine/2 - 1 or i >= ix_fine/2 + ix_low + 1: # for fine grid
            x[margin + ii] = x[margin + ii - 1] + dx_fine
        else: # for coarse grid
            x[margin + ii] = x[margin + ii - 1] + dx_low
        
    for i in range(0,margin): # for margin 
        x[i] = x[ixg - 2*margin + i] - (xmax - xmin)
        x[ixg - margin + i] = x[margin + i] + (xmax - xmin)

    
    return x

######################################################
######################################################
def upgrade_resolution(
        self,caseid,n
        ,xmin,xmax,ymin,ymax,zmin,zmax 
        ,ix,jx,kx
        ,endian='<',end_step=False
        ,memory_saving=False
        ,ix_ununi=32,dx00=4.8e6,x_ununif=False
        ):
    '''
    This function chabges the resolution and output for the next
    calculation.

    Parameters:
        self (R2D2_data): instance of R2D2_data
        caseid (str): caseid of destination directory e.g. 'd002'
        n (int): time step of upgrade data.
                 If end_step=True, n is ignored
        xmax (float): max location in x direction
        xmin (float): min location in x direction
        ymax (float): max location in y direction
        ymin (float): min location in y direction
        zmax (float): max location in z direction
        zmin (float): min location in z direction

        ix (int): updated grid point in x direction
        jx (int): updated grid point in y direction
        kx (int): updated grid point in z direction

        endian (str): endian
        end_step (bool): if true final time step is used for upgrade
                         and parameter n is ignored    
        memory_saving (bool): If true, upgraded variables are saved 
                              in memory separately for saving memory
        
        below prameters are effective only when x_ununif=True
        ix_ununi (int): number of grid in uniform grid region
        dx00 (float): grid spacing in uniform grid region
        x_ununif (bool): whethere ununiform grid is used

    '''
    import os
    import os.path
    import numpy as np
    import pyR2D2
    from . import common
    import sys
    from tqdm import tqdm

    self.up = {}

    ## read Model S based stratification
    ## data is stored in self.models dictionary
    self.models_init() 

    ## check if the origin and the destination are different
    if self.p['datadir'][:-6] == '../run/'+caseid:
        print('The origin and the destination are the same. stop..')
        return
    
    ## check if the destination directory
    ## if not this method finishes
    if os.path.exists('../run/'+caseid):
        print('The destination directory '+caseid+' already exist. stop..')
        return
    else:
        os.system('rsync -avP --exclude="data" '+self.p['datadir'][:-6]+'/'+' ../run/'+caseid)

    ## number of grid after upgrade
    self.up['ix'] = ix
    self.up['jx'] = jx
    self.up['kx'] = kx

    ## number of grid with margin after upgrade
    self.up['ixg'] = self.up['ix'] + 2*self.p['margin']
    self.up['jxg'] = self.up['jx'] + 2*self.p['margin']
    self.up['kxg'] = self.up['kx'] + 2*self.p['margin']
    
    if x_ununif:
        ## generate upgraded coordinate in ununiform geometry
        self.up['x'] = gen_coord_ununiform_top(xmax,xmin,self.up['ix'],self.p['margin'],dx00,ix_ununi)
    else:
        ## generate upgraded coordinate in uniform geometry
        self.up['x'] = gen_coord(xmax,xmin,self.up['ix'],self.p['margin'])

    ## generate upgraded coordinate
    self.up['y'] = gen_coord(ymax,ymin,self.up['jx'],self.p['margin'])
    self.up['z'] = gen_coord(zmax,zmin,self.up['kx'],self.p['margin'])
    
    ## read checkpoint data of original case
    self.read_qq_check(n,silent=True,end_step=end_step)

    ## prepare checkpoint data for upgrade data
    if memory_saving:
        self.qu = np.zeros((self.up['ixg'],self.up['jxg'],self.up['kxg']))
    else:
        self.qu = np.zeros((self.up['ixg'],self.up['jxg'],self.up['kxg'],self.p['mtype']))

    os.makedirs('../run/'+caseid+'/data/param/',exist_ok=True)
    os.makedirs('../run/'+caseid+'/data/qq/',exist_ok=True)
    os.makedirs('../run/'+caseid+'/data/slice/',exist_ok=True)
    os.makedirs('../run/'+caseid+'/data/remap/qq/',exist_ok=True)
    os.makedirs('../run/'+caseid+'/data/remap/vl/',exist_ok=True)
    os.makedirs('../run/'+caseid+'/data/time/mhd/',exist_ok=True)
    os.makedirs('../run/'+caseid+'/data/time/tau/',exist_ok=True)
    os.makedirs('../run/'+caseid+'/data/tau/',exist_ok=True)

    print('### Upgrade starts ###')
    for m in tqdm(range(0,self.p['mtype'])):
        if memory_saving:
            self.qu = pyR2D2.fortran_util.interp(self.p['xg'],self.p['yg'],self.p['zg'], \
                                                self.up['x'],self.up['y'],self.up['z'], \
                                                self.qc[:,:,:,m])
            self.qu.reshape([self.up['ixg']*self.up['jxg']*self.up['kxg']] \
                            ,order='F').astype(endian+'d').tofile('../run/'+caseid+'/data/qq/qq'+'{0:02d}'.format(m)+'.dac.e')            
        else:
            self.qu[:,:,:,m] = pyR2D2.fortran_util.interp(self.p['xg'],self.p['yg'],self.p['zg'], \
                                                  self.up['x'],self.up['y'],self.up['z'], \
                                                  self.qc[:,:,:,m])
            

    if not memory_saving:
        print('### Save data ###')
        self.qu.reshape([self.p['mtype']*self.up['ixg']*self.up['jxg']*self.up['kxg']] \
                        ,order='F').astype(endian+'d').tofile('../run/'+caseid+'/data/qq/qq.dac.e')


    def sign_judge(value):
        if np.sign(value) == 1.0:
            sign = "+"
        else:
            sign = "-"
        return sign

    def change_judge(value,self,key,enp=True):
        if value == self.p[key]:
            change = ' (unchange)'
        else:
            if enp:
                change = ' \033[31m(changed)\033[0m'
            else:
                change = ' (changed)'
                
        return change
    
    t = np.zeros(1)
    t.reshape([1] \
            ,order='F').astype(endian+'d').tofile('../run/'+caseid+'/data/time/mhd/t.dac.e')

    f = open('../run/'+caseid+'/data/param/nd.dac',mode='w')
    f.write(str(0).rjust(8)+str(0).rjust(8))
    f.close()

    f = open('../run/'+caseid+'/data/cont_log.txt',mode='w')
    f.write('The initial condition of this run is upgraded data from the other directory\n')
    f.write('This file describes the original and the upgraded data\n\n')
    f.write('### Original data ###\n\n')
    f.write('Server: '+self.p['server']+'\n')
    f.write('Datadir: '+self.p['datadir']+'\n')
    if end_step:
        cstep = str(self.p['nd'])
    else:
        cstep = str(n)
    f.write('Output step: '+cstep+'\n')

    value = self.p['xmin'] - self.p['rstar']
    f.write('xmin = rstar '+sign_judge(value)+'{:.4e}'.format(abs(value))+' or '
            +'{:.3f}'.format(self.p['xmin']/self.p['rstar'])+'*rstar'
            +change_judge(xmin,self,'xmin',enp=False)+'\n')
    value = self.p['xmax'] - self.p['rstar']
    f.write('xmax = rstar '+sign_judge(value)+'{:.4e}'.format(abs(value))+' or '
            +'{:.3f}'.format(self.p['xmax']/self.p['rstar'])+'*rstar'
            +change_judge(xmax,self,'xmax',enp=False)+'\n')
    if self.p['geometry'] == 'Spherical':
        print(self.p['ymin']/np.pi*180)
        f.write('ymin = '+'{:.4e}'.format(self.p['ymin']/np.pi*180)+' [rad] '+change_judge(ymin,self,'ymin',enp=False)+'\n')
        f.write('ymax = '+'{:.4e}'.format(self.p['ymax']/np.pi*180)+' [rad] '+change_judge(ymax,self,'ymax',enp=False)+'\n')
        f.write('zmin = '+'{:.4e}'.format(self.p['zmin']/np.pi*180)+' [rad] '+change_judge(zmin,self,'zmin',enp=False)+'\n')
        f.write('zmax = '+'{:.4e}'.format(self.p['zmax']/np.pi*180)+' [rad] '+change_judge(zmax,self,'zmax',enp=False)+'\n')
    else:
        f.write('ymin = '+'{:.4e}'.format(self.p['ymin'])+change_judge(ymin,self,'ymin',enp=False)+'\n')
        f.write('ymax = '+'{:.4e}'.format(self.p['ymax'])+change_judge(ymax,self,'ymax',enp=False)+'\n')
        f.write('zmin = '+'{:.4e}'.format(self.p['zmin'])+change_judge(zmin,self,'zmin',enp=False)+'\n')
        f.write('zmax = '+'{:.4e}'.format(self.p['zmax'])+change_judge(zmax,self,'zmax',enp=False)+'\n')
    f.write('\n')
    
    f.write('nx0*ix0 = '+str(self.p['ix'])+change_judge(self.up['ix'],self,'ix',enp=False)+'\n')
    f.write('ny0*jx0 = '+str(self.p['jx'])+change_judge(self.up['jx'],self,'jx',enp=False)+'\n')
    f.write('nz0*kx0 = '+str(self.p['kx'])+change_judge(self.up['kx'],self,'kx',enp=False)+'\n')

    if x_ununif:
        uniform_flag = '.true.'
    else:
        uniform_flag = '.false.'

    f.write('uniform_flag = '+str(self.p['ununiform_flag']) \
            +change_judge(x_ununif,self,'ununiform_flag',enp=False)+'\n')

    f.write('ix_ununi = '+str(self.p['ix_ununi'])+change_judge(ix_ununi,self,'ix_ununi',enp=False)+'\n')
    f.write('dx00 = ''{:.4e}'.format(self.p['dx00'])+change_judge(dx00,self,'dx00',enp=False)+'\n')

    f.write('\n')
    f.write('###############################################')
    f.write('\n')
    f.write('###############################################')
    f.write('\n')
    f.write('###############################################')
    f.write('\n')
    f.write('### Upgrade data ###')
    f.write('\n\n')

    value = xmin - self.p['rstar']
    f.write('xmin = rstar '+sign_judge(value)+'{:.4e}'.format(abs(value))+' or '
            +'{:.3f}'.format(xmin/self.p['rstar'])+'*rstar'
            +change_judge(xmin,self,'xmin',enp=False)+'\n')
    value = xmax - self.p['rstar']
    f.write('xmax = rstar '+sign_judge(value)+'{:.4e}'.format(abs(value))+' or '
            +'{:.3f}'.format(xmax/self.p['rstar'])+'*rstar'
            +change_judge(xmax,self,'xmax',enp=False)+'\n')

    if self.p['geometry'] == 'Spherical':
        f.write('ymin = '+'{:.4e}'.format(ymin/np.pi*180)+' [rad] '+change_judge(ymin,self,'ymin',enp=False)+'\n')
        f.write('ymax = '+'{:.4e}'.format(ymax/np.pi*180)+' [rad] '+change_judge(ymax,self,'ymax',enp=False)+'\n')
        f.write('zmin = '+'{:.4e}'.format(zmin/np.pi*180)+' [rad] '+change_judge(zmin,self,'zmin',enp=False)+'\n')
        f.write('zmax = '+'{:.4e}'.format(zmax/np.pi*180)+' [rad] '+change_judge(zmax,self,'zmax',enp=False)+'\n')
    else:
        f.write('ymin = '+'{:.4e}'.format(ymin)+change_judge(ymin,self,'ymin',enp=False)+'\n')
        f.write('ymax = '+'{:.4e}'.format(ymax)+change_judge(ymax,self,'ymax',enp=False)+'\n')
        f.write('zmin = '+'{:.4e}'.format(zmin)+change_judge(zmin,self,'zmin',enp=False)+'\n')
        f.write('zmax = '+'{:.4e}'.format(zmax)+change_judge(zmax,self,'zmax',enp=False)+'\n')
        
        f.write('\n')
    
    f.write('nx0*ix0 = '+str(self.up['ix'])+change_judge(self.up['ix'],self,'ix',enp=False)+'\n')
    f.write('ny0*jx0 = '+str(self.up['jx'])+change_judge(self.up['jx'],self,'jx',enp=False)+'\n')
    f.write('nz0*kx0 = '+str(self.up['kx'])+change_judge(self.up['kx'],self,'kx',enp=False)+'\n')

    if x_ununif:
        uniform_flag = '.true.'
    else:
        uniform_flag = '.false.'

    f.write('uniform_flag = '+str(x_ununif) \
            +change_judge(x_ununif,self,'ununiform_flag',enp=False)+'\n')

    f.write('ix_ununi = '+str(ix_ununi)+change_judge(ix_ununi,self,'ix_ununi',enp=False)+'\n')
    f.write('dx00 = ''{:.4e}'.format(dx00)+change_judge(dx00,self,'dx00',enp=False)+'\n')

    
    f.close()
    
    print(' ')
    print('### Data upgrade funished ###')
    print('Please use following parameters')
    print(' ')
    print('caseid is \033[31m'+caseid+'\033[0m')
    print(' ')
        
    value = xmin - self.p['rstar']
    print('xmin = rstar '+sign_judge(value),'{:.4e}'.format(abs(value)),' or '
          ,'{:.3f}'.format(xmin/self.p['rstar'])+'*rstar'
          ,change_judge(xmin,self,'xmin'))
    value = xmax - self.p['rstar']
    print('xmax = rstar '+sign_judge(value),'{:.4e}'.format(abs(value)),' or '
          ,'{:.3f}'.format(xmax/self.p['rstar'])+'*rstar'
          ,change_judge(xmax,self,'xmax'))

    if self.p['geometry'] == 'Spherical':
        print('ymin = '+'{:.4e}'.format(ymin/np.pi*180),change_judge(ymin,self,'ymin'))
        print('ymax = '+'{:.4e}'.format(ymax/np.pi*180),change_judge(ymax,self,'ymax'))
        print('zmin = '+'{:.4e}'.format(zmin/np.pi*180),change_judge(zmin,self,'zmin'))
        print('zmax = '+'{:.4e}'.format(zmax/np.pi*180),change_judge(zmax,self,'zmax'))
    else:
        print('ymin = '+'{:.4e}'.format(ymin),change_judge(ymin,self,'ymin'))
        print('ymax = '+'{:.4e}'.format(ymax),change_judge(ymax,self,'ymax'))
        print('zmin = '+'{:.4e}'.format(zmin),change_judge(zmin,self,'zmin'))
        print('zmax = '+'{:.4e}'.format(zmax),change_judge(zmax,self,'zmax'))

    print(' ')
    
    print('nx0*ix0 = '+str(self.up['ix']),change_judge(self.up['ix'],self,'ix'))
    print('ny0*jx0 = '+str(self.up['jx']),change_judge(self.up['jx'],self,'jx'))
    print('nz0*kx0 = '+str(self.up['kx']),change_judge(self.up['kx'],self,'kx'))

    print(' ')

    if x_ununif:
        uniform_flag = '.true.'
    else:
        uniform_flag = '.false.'

    print('uniform_flag = '+uniform_flag \
          ,change_judge(x_ununif,self,'ununiform_flag'))

    if x_ununif:
        print('ix_ununi = '+str(ix_ununi),change_judge(ix_ununi,self,'ix_ununi'))
        print('dx00 = ''{:.4e}'.format(dx00),change_judge(dx00,self,'dx00'))
    else:
        print('You plan to use uniform grid spacing.')
        print('Do not care about ix_ununi and dx00 in geometry_def.F90')
