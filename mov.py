import numpy as np
from tqdm import tqdm
import R2D2
import sys, os
import matplotlib.pyplot as plt
import mov_util

caseid = R2D2.util.caseid_select(locals())
datadir="../run/"+caseid+"/data/"

R2D2.util.initialize_instance(locals(),'d')
d = R2D2.R2D2_data(datadir,self_old=d)
R2D2.util.locals_define(d,locals())

pngdir="../figs/"+caseid+"/mov/"
os.makedirs(pngdir,exist_ok=True)

n0 = R2D2.util.define_n0(d,locals(),nd_type='nd')

print("Maximum time step= ",d.p['nd']," time ="\
          ,d.p['dtout']*float(d.p['nd'])/3600./24.," [day]")

plt.close('all')

if not d.p['geometry'] == 'Cartesian':
    d.p['zz'],d.p['yy'] = np.meshgrid(d.p['z'],d.p['y']-0.5*np.pi)
    
    xe = np.zeros(d.p['ix']+1)
    ye = np.zeros(d.p['jx']+1)

    xe[0] = d.p['xmin']
    xe[d.p['ix']] = d.p['xmax']
    ye[0] = 0.e0
    ye[d.p['jx']] = np.pi
    for i in range(1,d.p['ix']):
        xe[i] = xe[i-1] + 2*(d.p['x'][i-1] - xe[i-1])

    for j in range(1,d.p['jx']):
        ye[j] = ye[j-1] + 2*(d.p['y'][j-1] - ye[j-1])

    RAE, THE = np.meshgrid(xe,ye,indexing='ij')
    d.p['XX'], d.p['YY'] = RAE*np.cos(THE), RAE*np.sin(THE)

# read initial time
t0 = d.read_time(0,silent=True)

for n in tqdm(range(n0,d.p['nd']+1)):
    # read data
    t = d.read_time(n,silent=True) 
    d.read_vc(n,silent=True)
    
    ##############################
    
    if n == n0:
        tight_layout_flag = True
    else:
        tight_layout_flag = False
    
    if d.p['geometry'] == 'Cartesian':
        d.read_qq_tau(n*int(d.p['ifac']),silent=True)
        tu_height = d.qt['he'][d.p['jc'],:]
        
        inm = d.qt['in'].mean() # mean intensity
        inrms = np.sqrt(((d.qt['in'] - inm)**2).mean()) # RMS intensity
        frms = 2.

        sem, tmp = np.meshgrid(d.vc['sem'].mean(axis=1),z,indexing='ij')
        serms, tmp = np.meshgrid(np.sqrt((d.vc['serms']**2).mean(axis=1)),z,indexing='ij')
        
        bb = np.sqrt(d.vc["bx_xz"]**2 + d.vc["by_xz"]**2 + d.vc["bz_xz"]**2)        

        vls = [d.qt['in']*1.e-10,
               d.qt['bx'],
               (d.vc['se_xz'] - sem)/serms,
               bb
               ]
        vmaxs = [(inm+inrms*frms)*1.e-10, # intensity
                 2.5e3, # LoS B
                 2., # entropy
                 8.e3, # magnetic field strength
                 ]
        vmins = [(inm-inrms*frms)*1.e-10, # intensity
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
        d.read_qq_select(d.p['xmax'],n,silent=True)
        vxrms = np.sqrt((d.qs['vx']**2).mean())
        bxrms = np.sqrt((d.qs['bx']**2).mean())
        
        sem, tmp   = np.meshgrid((d.vc['sem']*SINY).sum(axis=1)/SINYM,d.p['y'],indexing='ij')
        serms, tmp = np.meshgrid(np.sqrt((d.vc['serms']**2*SINY).sum(axis=1)/SINYM),d.p['y'],indexing='ij')
        
        if serms.max() != 0:
            se_value = (d.vc['se_xy']-sem)/serms
        else:
            se_value = np.zeros((d.p['ix'],d.p['jx']))
        
        vls = [d.qs['vx']*1.e-2,
               d.qs['bx']*1.e-3,
               se_value,
               d.vc['bzm']
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
    if(n != d.p['nd']):
        plt.clf()
        

        
    
 
