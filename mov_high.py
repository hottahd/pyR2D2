import os

import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm

import pyR2D2
import mov_util

caseid = pyR2D2.util.caseid_select(locals())
datadir="../run/"+caseid+"/data/"

pyR2D2.util.initialize_instance(locals(),'d')
d = pyR2D2.Data(datadir,self_old=d)
pyR2D2.util.locals_define(d,locals())

pngdir="../figs/"+caseid+"/mov_high/"
os.makedirs(pngdir,exist_ok=True)

n0 = pyR2D2.util.define_n0(d,locals(),nd_type='nd_tau')

print("Maximum time step= ",d.nd_tau," time ="\
          ,d.dtout_tau*float(d.nd_tau)/3600./24.," [day]")

plt.close('all')

if d.geometry == 'Cartesian':
    TE0, tmp = np.meshgrid(d.te0, d.z, indexing='ij')
else:
    d.zz, d.yy = np.meshgrid(d.z, d.y - 0.5*np.pi)
    
    xe = np.zeros(d.ix + 1)
    ye = np.zeros(d.jx + 1)

    xe[0] = d.xmin
    xe[d.ix] = d.xmax
    ye[0] = 0.e0
    ye[d.jx] = np.pi
    for i in range(1,d.ix):
        xe[i] = xe[i-1] + 2*(d.x[i-1] - xe[i-1])

    for j in range(1,d.p['jx']):
        ye[j] = ye[j-1] + 2*(d.y[j-1] - ye[j-1])

    RAE, THE = np.meshgrid(xe,ye,indexing='ij')
    d.p.XX, d.p.YY = RAE*np.cos(THE), RAE*np.sin(THE)

# read initial time
t0 = d.time_read(0,verbose=False)

for n in tqdm(range(n0,d.nd_tau + 1)):
    # read data
    t = d.time_read(n,verbose=False,tau=True)
    
    ##############################
    
    if n == n0:
        tight_layout_flag = True
    else:
        tight_layout_flag = False
    
    if d.geometry == 'Cartesian':
        d.qt.read(n)
        d.qs.read(np.argmin(abs(d.y_slice - 0.5*d.ymax)),'y',n)
        tu_height = d.qt.he[np.argmax(d.y > d.y_slice[np.argmin(abs(d.y_slice - 0.5*d.ymax))]),:]
                
        rtm = d.qt.rt.mean() # mean intensity
        rtrms = np.sqrt(((d.qt.rt - rtm)**2).mean()) # RMS intensity
        frms = 2.
        
        bb = np.sqrt(d.qs.bx**2 + d.qs.by**2 + d.qs.bz**2)

        vls = [d.qt.rt*1.e-10,
               d.qt.bx,
               d.qs.te + TE0,
               bb
               ]
        vmaxs = [(rtm + rtrms*frms)*1.e-10, # intensity
                 2.5e3, # LoS B
                 d.te0.max(),
                 8.e3, # magnetic field strength
                 ]
        vmins = [(rtm - rtrms*frms)*1.e-10, # intensity
                 -2.5e3, # LoS B
                 d.te0.min(),
                 0 # magnetic field strength
                 ]
        titles = ['Emergent intensity\n'+r' $\mathrm{[10^{10}~erg~cm^{-2}~ster^{-1}~s^{-1}]}$',
                  "LOS magnetic field\n"+r"at $\tau=1~\mathrm{[G]}$",
                  r"$T~\mathrm{[K]}$",
                  r"$|B|~\mathrm{[G]}$"
                ]

        mov_util.mov_cartesian_photo_2x2(d,t-t0,vls,tu_height,vmaxs,vmins,titles,tight_layout_flag=tight_layout_flag)
    else: # Spherical geometry including Yin-Yang
        d.qs.read(np.argmin(abs(d.x_slice - d.xmax)), 'x', n, silent=True) # xmaxに一番近いところ
        vxrms = np.sqrt((d.qs.vx_yin**2).mean())
        bxrms = max(np.sqrt((d.qs.bx_yin**2).mean()),1e-2)
    
        vfac = 1.e-2
        bfac = 1.e-3
    
        vls = [
            {'Yin': d.qs.vx_yin*vfac, 'Yan': d.qs.vx_yan*vfac},
            {'Yin': d.qs.bx_yin*bfac, 'Yan': d.qs.bx_yan*bfac},
        ]
        
        vmaxs = [
            2*vxrms*vfac,
            2*bxrms*bfac,
        ]
        
        vmins = [
            -vmaxs[0],
            -vmaxs[1]
        ]
        
        titles = [r'Radial velocity $v_r~\mathrm{[m~s^{-1}]}$', r'Radial magnetic field $B_r~\mathrm{[kG]}$']
        mov_util.mov_yinyang_2(d,t,vls,vmaxs,vmins,titles,tight_layout_flag=tight_layout_flag)
        
    plt.pause(0.1)
    plt.savefig(pngdir+"py"+'{0:08d}'.format(n)+".png")
    if(n != d.nd_tau):
        plt.clf()
        

        
    
 
