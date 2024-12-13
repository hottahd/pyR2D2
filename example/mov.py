import os

import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import mov_util
import pyR2D2

caseid = pyR2D2.util.caseid_select(locals())
datadir="../run/"+caseid+"/data/"

pyR2D2.util.initialize_instance(locals(),'d')
d = pyR2D2.Data(datadir,self_old=d)
pyR2D2.util.locals_define(d,locals())

pngdir="../figs/"+caseid+"/mov/"
os.makedirs(pngdir,exist_ok=True)

n0 = pyR2D2.util.define_n0(d,locals(),nd_type='nd')

print("Maximum time step= ",d.nd," time ="\
          ,d.dtout*float(d.nd)/3600./24.," [day]")

plt.close('all')

if not d.geometry == 'Cartesian':
    d.zz, d.yy = np.meshgrid(d.z, d.y - 0.5*np.pi)
    
    xe = np.zeros(d.ix + 1)
    ye = np.zeros(d.jx + 1)

    xe[0] = d.xmin
    xe[d.ix] = d.xmax
    ye[0] = 0.e0
    ye[d.jx] = np.pi
    for i in range(1,d.ix):
        xe[i] = xe[i-1] + 2*(d.x[i-1] - xe[i-1])

    for j in range(1,d.jx):
        ye[j] = ye[j-1] + 2*(d.y[j-1] - ye[j-1])

    RAE, THE = np.meshgrid(xe,ye,indexing='ij')
    d.p.XX, d.p.YY = RAE*np.cos(THE), RAE*np.sin(THE)

    X, Y = np.meshgrid(d.x,d.y,indexing='ij')
    SINY = np.sin(Y)
    SINYM = SINY.sum(axis=1)

# read initial time
t0 = d.time_read(0,verbose=False)

for n in tqdm(range(n0,d.nd + 1)):
    # read data
    t = d.time_read(n,verbose=False)
    d.vc.read(n)
    
    ##############################
    
    if n == n0:
        tight_layout_flag = True
    else:
        tight_layout_flag = False
    
    if d.geometry == 'Cartesian':
        d.qt.read(n*int(d.ifac))
        tu_height = d.qt.he[d.jc, :]
        
        rtm = d.qt.rt.mean() # mean intensity
        rtrms = np.sqrt(((d.qt.rt - rtm)**2).mean()) # RMS intensity
        frms = 2.

        sem, tmp = np.meshgrid(d.vc.sem.mean(axis=1),z,indexing='ij')
        serms, tmp = np.meshgrid(np.sqrt((d.vc.serms**2).mean(axis=1)),d.z,indexing='ij')
        
        bb = np.sqrt(d.vc.bx_xz**2 + d.vc.by_xz**2 + d.vc.bz_xz**2)

        vls = [d.qt.rt*1.e-10,
               d.qt.bx,
               (d.vc.se_xz - sem)/serms,
               bb
               ]
        vmaxs = [(rtm + rtrms*frms)*1.e-10, # intensity
                 2.5e3, # LoS B
                 2., # entropy
                 8.e3, # magnetic field strength
                 ]
        vmins = [(rtm - rtrms*frms)*1.e-10, # intensity
                 -2.5e3, # LoS B
                 -2., # entropy
                 0 # magnetic field strength
                 ]
        titles = ['Emergent intensity\n'+r' $\mathrm{[10^{10}~erg~cm^{-2}~ster^{-1}~s^{-1}]}$',
                  "LOS magnetic field\n"+r"at $\tau=1~\mathrm{[G]}$",
                  r"$\left(s-\langle s\rangle\right)/s_\mathrm{RMS}$",
                  r"$|B|~\mathrm{[G]}$"
                ]
        mov_util.mov_cartesian_photo_2x2(d,t,vls,tu_height,vmaxs,vmins,titles,tight_layout_flag=tight_layout_flag)
    else: # Spherical geometry including Yin-Yang
        d.qx.read(d.xmax, n)
        vxrms = np.sqrt((d.qx.vx**2).mean())
        bxrms = np.sqrt((d.qx.bx**2).mean())
        
        sem, tmp   = np.meshgrid((d.vc.sem*SINY).sum(axis=1)/SINYM,d.y,indexing='ij')
        serms, tmp = np.meshgrid(np.sqrt((d.vc.serms**2*SINY).sum(axis=1)/SINYM),d.y,indexing='ij')
        
        if serms.max() != 0:
            se_value = (d.vc.se_xy - sem)/serms
        else:
            se_value = np.zeros((d.ix, d.jx))
        
        vls = [d.qx.vx*1.e-2,
               d.qx.bx*1.e-3,
               se_value,
               d.vc.bzm
               ]
        vmaxs = [2*vxrms*1.e-2, # radial velocity
                 2*bxrms*1.e-3, # radial magnetic field
                 2., # entropy
                 8000, # magnetic field strength
                 ]
        vmins = [-vmaxs[0], # radial velocity
                 -vmaxs[1], # radial magnetic field
                 -2., # entropy
                 -8000 # magnetic field strength
                 ]
        titles = [r'Radial velocity $v_r~[\mathrm{m~s^{-1}}]$',
                  r'Radial magnetic field $B_r~[\mathrm{kG}]$',
                  r"$\left(s-\langle s\rangle\right)/s_\mathrm{RMS}$",
                  r"$\langle B_\phi\rangle~\mathrm{[G]}$"
                ]
                
        mov_util.mov_spherical_2x2(d,t,vls,vmaxs,vmins,titles,tight_layout_flag=tight_layout_flag)
    plt.pause(0.1)
    plt.savefig(pngdir+"py"+'{0:08d}'.format(n)+".png")
    if(n != d.nd):
        plt.clf()
        

        
    
 
