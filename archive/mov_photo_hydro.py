import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from tqdm import tqdm
import pyR2D2
import sys
import os

try:
    caseid
except NameError:
    print("input caseid id (3 digit)")
    caseid = 0
    caseid = input()
    caseid = "d"+caseid.zfill(3)
    
datadir="../run/"+caseid+"/data/"
pngdir="../figs/"+caseid+"/mov_photo_hydro/"
os.makedirs(pngdir,exist_ok=True)

d = pyR2D2.R2D2_data(datadir)
for key in d.p:
    exec('%s = %s%s%s' % (key, 'd.p["',key,'"]'))
    
try:
    n0
except NameError:
    n0 = 0
if  n0 > d.p["nd_tau"]:
    n0 = d.p["nd_tau"]

print("Maximum time step (nd_tau) = ",nd_tau," time ="\
          ,dtout*float(nd)/3600./24.," [day]")

plt.close('all')

# read initial time
t0 = d.read_time(0,silent=True)

yran = ymax - ymin
xran = min(xmax-xmin,yran)

xsize = 6
ysize = xsize*(yran + xran)/yran
fig = plt.figure(num='mov_photo_hydro',figsize=(xsize + 0.5,ysize))

grid = GridSpec(2,1,height_ratios=[yran,xran])

te2, tmp = np.meshgrid(te0,y,indexing='ij')

#n0 = 18
#nd_tau = n0

for n in tqdm(range(n0,nd_tau+1)):
#for n in range(0,1):
    #print(n)
    ##############################
    # read time
    t = d.read_time(n,tau=True,silent=True)
        
    ##############################
    # read time
    d.read_qq_tau(n,silent=True)

    ##############################
    # read value

    k_slice = 1
    d.read_qq_slice(k_slice,'z',n,silent=True)
    kl = np.argmin(np.abs(z - z_slice[k_slice]))
    ##############################

    shading = "auto"

    lfac = 1.e-8
    
    ax1 = fig.add_subplot(grid[0,0],aspect='equal')
    ax2 = fig.add_subplot(grid[1,0],aspect='equal')
    
    in_qs = 2.25e10
    in0 = d.qt["in"].copy()/in_qs
    in0s = np.roll(in0,[jx//2-jc,kx//2-kc],axis=[0,1])
    im1 = ax1.pcolormesh(y*lfac,z*lfac,in0s.transpose(),cmap='gist_gray',vmax=1.4,vmin=0.5,shading=shading)
    ax1.set_title(r"Emergent intensity $I/I_\odot$")
    ax2.set_title(r"$T$ [$10^3$ K]")
    im2 = ax2.pcolormesh(y*lfac,(x-rsun)*lfac,(d.ql['te']+te2)*1.e-3,vmin=2.,vmax=15,cmap='rainbow',shading=shading)
    ax2.plot(y*lfac,(d.qt['he'][:,kl]-rsun)*lfac,color='w')

    fig.subplots_adjust(hspace=0.1,wspace=0.05,top=0.96,bottom=0.1,left=0.1,right=0.9)
    box1 = ax1.get_position().bounds
    box2 = ax2.get_position().bounds
    cax_left = 0.905
    cax_width = 0.02
    cax1 = fig.add_axes([cax_left,box1[1],cax_width,box1[3]])
    cax2 = fig.add_axes([cax_left,box2[1],cax_width,box2[3]])

    cl1 = fig.colorbar(im1,cax=cax1)
    cl2 = fig.colorbar(im2,cax=cax2)
    
    ax1.tick_params(labelbottom=False)
    ax2.set_xlabel('Mm')
    for ax in [ax1,ax2]:
        ax.set_ylabel('Mm')

    ax2.annotate(text="$t$="+"{:.2f}".format((t-t0)/60/60)+" [hour]"\
                     ,xy=[0.04,0.02],xycoords="figure fraction"\
                     ,color='black',fontsize=24)#,bbox=bbox_props)

#    if(n == n0):
#        fig.tight_layout(pad=0.1)
    
    plt.pause(0.1)
    plt.savefig(pngdir+"py"+'{0:08d}'.format(n)+".png")

    if(n != nd_tau):
        plt.clf()
