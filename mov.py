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

yran = ymax - ymin
xran = min(xmax-xmin,yran)

xsize = 12
ysize = xsize*(yran + xran)/2/yran
fig = plt.figure(num=1,figsize=(xsize,ysize))

grid = GridSpec(2,2,height_ratios=[yran,xran])

for n in tqdm(range(n0,nd+1)):
    ##############################
    # read time
    t = d.read_time(n,silent=True)
        
    ##############################
    # read time
#    if xmax > rsun:
    d.read_qq_tau(n*int(ifac),silent=True)

    ##############################
    # read value
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
    ax1.pcolormesh(y*lfac,z*lfac,d.qt['in'].T,cmap='inferno',vmax=inm+inrms*frms,vmin=inm-inrms*frms,shading=shading)
    ax2.pcolormesh(y*lfac,z*lfac,d.qt['bx'].T,cmap='gist_gray',vmax=2.5e3,vmin=-2.5e3,shading=shading)
    ax3.pcolormesh(y*lfac,(x-rstar)*lfac,(d.vc['se_xy'] - d.vc['sem'])/d.vc['serms'],vmin=-2.,vmax=2.,shading=shading)
    ax3.contour(y*lfac,(x-rstar)*lfac,d.vc['tu_xy'],levels=[1.],colors="w")
    
    bb = np.sqrt(d.vc["bx_xy"]**2 + d.vc["by_xy"]**2 + d.vc["bz_xy"]**2)
    ax4.pcolormesh(y*lfac,(x-rstar)*lfac,bb,vmax=8.e3,vmin=0.,cmap='gist_heat',shading=shading)
    
    ax4.contour(y*lfac,(x-rstar)*lfac,d.vc['tu_xy'],levels=[1.],colors="w")

    for ax in [ax1,ax2]:
        ax.tick_params(labelbottom=False)
    for ax in [ax2,ax4]:
        ax.tick_params(labelleft=False)
    for ax in [ax3,ax4]:
        ax.set_xlabel('$y$ [Mm]')
    for ax in [ax1,ax3]:
        ax.set_ylabel('$x$ [Mm]')
    
    for ax in [ax3,ax4]:
        ax.set_ylim((max(xmax-yran,xmin)-rstar)*lfac,(xmax-rstar)*lfac)
        
    ax1.set_title('Emergent intensity')        
    ax2.set_title(r"LOS magnetic field at $\tau=1$")
    ax3.set_title(r"$(s-\langle s\rangle)/s_\mathrm{RMS}$")
    ax4.set_title(r"$|B|$")
    

    ax3.annotate(text="$t$="+"{:.2f}".format((t-t0)/60/60)+" [hour]"\
                     ,xy=[0.03,0.03],xycoords="figure fraction"\
                     ,color='black')#,bbox=bbox_props)

    if(n == n0):
        fig.tight_layout(pad=0.1)
    
    plt.pause(0.1)
    plt.savefig(pngdir+"py"+'{0:08d}'.format(n)+".png")

    if(n != nd):
        plt.clf()
