import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patheffects import withStroke
from tqdm import tqdm
import R2D2
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
pngdir="../figs/"+caseid+"/mov_photo/"
os.makedirs(pngdir,exist_ok=True)

d = R2D2.R2D2_data(datadir)
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

xsize = 16
ysize = 8
fig = plt.figure('mov_photo',figsize=(xsize,ysize))

n0_shift = 0
#n0 = 1720
#nd_tau = n0

d.read_vc(0,silent=True)

for n in tqdm(range(n0,nd_tau+1)):
    #print(n)
    ##############################
    # read time
    t = d.read_time(n,tau=True,silent=True)
        
    ##############################
    # read time
    d.read_qq_tau(n,silent=True)

    ##############################

    shading = "auto"

    lfac = 1.e-8    
    ax1 = fig.add_subplot(121,aspect='equal')
    ax2 = fig.add_subplot(122,aspect='equal')
    
    in0 = d.qt["in"].copy()
    in0s = np.roll(in0,[jx//2-jc,kx//2-kc],axis=[0,1])
    ax1.pcolormesh(y*lfac,z*lfac,in0s.transpose(),cmap='gist_gray',vmax=3.3e10,vmin=1.4e10,shading=shading)
    ax1.set_ylabel("Mm")
    ax1.set_title("Emergent intensity")

    bx = np.roll(d.qt["bx"],[jx//2-jc,kx//2-kc],axis=[0,1])
    ax2.tick_params(labelleft=False)
    ax2.pcolormesh(y*lfac,z*lfac,bx.transpose(),cmap='gist_gray',vmax=3.e3,vmin=-3.e3,shading=shading)
    ax2.set_title(r"LOS magnetic field@$\tau=1$")

    for ax in [ax1,ax2]:
        ax.set_xlabel('Mm')

    if(n == n0):
        plt.tight_layout()

    ax1.annotate(text="$t$="+"{:.2f}".format((t-t0)/60/60)+" [hour]"\
                     ,xy=[5,5],xycoords="data",fontsize=25 \
                     ,color='white',path_effects=[withStroke(foreground='black',linewidth=3)])
        
    plt.pause(0.1)
    plt.savefig(pngdir+"py"+'{0:08d}'.format(n-n0_shift)+".png")

    if(n != nd_tau):
        plt.clf()
