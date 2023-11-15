import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from tqdm import tqdm
import R2D2
import sys
import os

caseid = R2D2.util.caseid_select(locals())
datadir="../run/"+caseid+"/data/"

R2D2.util.initialize_instance(locals(),'d')
d = R2D2.R2D2_data(datadir,self_old=d)
R2D2.util.locals_define(d,locals())

pngdir="../figs/"+caseid+"/mov/"
os.makedirs(pngdir,exist_ok=True)

R2D2.util.define_n0(d,locals())

print("Maximum time step= ",nd," time ="\
          ,dtout*float(nd)/3600./24.," [day]")

plt.close('all')

# read initial time
t0 = d.read_time(0,silent=True)

zran = zmax - zmin
yran = ymax - ymin
xran = min(xmax-xmin,zran)

xsize = 18
ysize = xsize*(yran + xran)/2/zran
ysize = ysize + min(xsize/ysize/1.5,2)
ysize_limit = 11
if ysize > ysize_limit:
    xsize = xsize/ysize*ysize_limit
    ysize = ysize_limit
fig = plt.figure(num=1,figsize=(xsize,ysize))

grid = GridSpec(2,2,height_ratios=[yran,xran])

for n in tqdm(range(n0,nd+1)):
    # read data
    t = d.read_time(n,silent=True) 
    d.read_qq_tau(n*int(ifac),silent=True)
    d.read_vc(n,silent=True)
    
    ##############################

    shading = "auto"

    lfac = 1.e-8
    
    ax1 = fig.add_subplot(grid[0,0],aspect='equal')
    ax2 = fig.add_subplot(grid[0,1],aspect='equal')
    ax3 = fig.add_subplot(grid[1,0],aspect='equal')
    ax4 = fig.add_subplot(grid[1,1],aspect='equal')
    
    inm = d.qt['in'].mean()
    inrms = np.sqrt(((d.qt['in'] - inm)**2).mean())
    frms = 2.
    
    # pcolormesh    
    im1 = ax1.pcolormesh(z*lfac,y*lfac,d.qt['in']*1.e-10,cmap='inferno',vmax=(inm+inrms*frms)*1.e-10,vmin=(inm-inrms*frms)*1.e-10,shading=shading)
    im2 = ax2.pcolormesh(z*lfac,y*lfac,d.qt['bx'],cmap='gist_gray',vmax=2.5e3,vmin=-2.5e3,shading=shading)
    sem, tmp = np.meshgrid(d.vc['sem'].mean(axis=1),z,indexing='ij')
    serms, tmp = np.meshgrid(np.sqrt((d.vc['serms']**2).mean(axis=1)),z,indexing='ij')
    
    im3 = ax3.pcolormesh(z*lfac,(x-rstar)*lfac,(d.vc['se_xz'] - sem)/serms,vmin=-2.,vmax=2.,shading=shading)
    ax3.contour(z*lfac,(x-rstar)*lfac,d.vc['tu_xz'],levels=[1.],colors="w")
    
    bb = np.sqrt(d.vc["bx_xz"]**2 + d.vc["by_xz"]**2 + d.vc["bz_xz"]**2)
    im4 = ax4.pcolormesh(z*lfac,(x-rstar)*lfac,bb,vmax=8.e3,vmin=0.,cmap='gist_heat',shading=shading)
    
    ax4.contour(z*lfac,(x-rstar)*lfac,d.vc['tu_xz'],levels=[1.],colors="w")

    for ax in [ax1,ax2]:
        ax.tick_params(labelbottom=False)
    for ax in [ax2,ax4]:
        ax.tick_params(labelleft=False)
    for ax in [ax3,ax4]:
        ax.set_xlabel('$z$ [Mm]')
    ax1.set_ylabel('$y$ [Mm]')
    ax3.set_ylabel('$x$ [Mm]')
    
    for ax in [ax3,ax4]:
        ax.set_ylim((max(xmax-yran,xmin)-rstar)*lfac,(xmax-rstar)*lfac)
    
    titles = ['Emergent intensity'+r' $\mathrm{[10^{10}~erg~cm^{-2}~ster^{-1}~s^{-1}]}$',
              r"LOS magnetic field at $\tau=1~\mathrm{[G]}$",
              r"$\left(s-\langle s\rangle\right)/s_\mathrm{RMS}$",
              r"$|B|~\mathrm{[G]}$"
    ]
    for ax, title in zip([ax1,ax2,ax3,ax4],titles):
        title = ax.set_title(title)
        title.set_position([0.01,1.02])
        title.set_ha('left')
        
    if(n == n0):
        fig.tight_layout()
        
    ax3.annotate(text="$t="+"{:.2f}".format((t-t0)/60/60)+"~\mathrm{[hour]}$"\
                     ,xy=[0.01,0.02],xycoords="figure fraction"\
                     ,color='black',fontsize=20)
    
    # add color bar
    box3 = ax3.get_position().bounds
    wid3, hei3 = box3[2], box3[3]
    
    color_bar_width = 0.12
    color_bar_height = 0.08/ysize
    
    for ax, im in zip([ax1,ax2,ax3,ax4],[im1,im2,im3,im4]):
        box = ax.get_position().bounds
        cax = fig.add_axes([box[0]+ box[2] - color_bar_width
                        ,box[1] + box[3] + color_bar_height*4
                        ,color_bar_width
                        ,color_bar_height
                        ])
        cbar = fig.colorbar(im,cax=cax,orientation='horizontal')
        cbar.ax.tick_params(labelsize=12)
    
    
    plt.pause(0.1)
    plt.savefig(pngdir+"py"+'{0:08d}'.format(n)+".png")

    if(n != nd):
        plt.clf()
