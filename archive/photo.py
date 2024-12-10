import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
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

dir="../run/"+caseid+"/data/"
pngdir="../figs/"+caseid+"/photo/"
os.makedirs(pngdir,exist_ok=True)

d = pyR2D2.R2D2_data(dir)
for key in d.p:
    exec('%s = %s%s%s' % (key, 'd.p["',key,'"]'))

try:
    n0
except NameError:
    n0 = 0
if  n0 > d.p["nd_tau"]:
    n0 = d.p["nd_tau"]

print("Maximum time step= ",nd_tau," time ="\
      ,dtout/ifac*float(nd_tau)/3600./24.," [day]")

plt.close('all')

xsize = 21
ysize = 7
fig = plt.figure(100,figsize=(xsize,ysize))

t0 = d.read_time(0,tau=True,silent=True)

#n0 = 0
#nd_tau = 500

#n0 = 801
#nd_tau = n0

for n in range(n0,nd_tau+1):
#for n in range(20,21):
    print(n)
    ##############################
    # read time
    t = d.read_time(n,tau=True,silent=True)

    ##############################
    # read time
    d.read_qq_tau(n,silent=True)

    shading = "auto"
    #shading = "groroud"
    ax1 = fig.add_subplot(131,aspect="equal")
    ax2 = fig.add_subplot(132,aspect="equal")
    ax3 = fig.add_subplot(133,aspect="equal")
    
    in0 = np.roll(d.qt["in"],[jx//2-jc,kx//2-kc],axis=(0,1))
    bx0 = np.roll(d.qt["bx001"],[jx//2-jc,kx//2-kc],axis=(0,1))
    bh0 = np.roll(sqrt(d.qt['by001']**2 + d.qt['bz001']**2),[jx//2-jc,kx//2-kc],axis=(0,1))

    lfac = 1.e-8
    ax1.pcolormesh(y*lfac,z*lfac,in0.transpose(),cmap='gist_gray',vmax=3.0e10,vmin=0.2e10,shading=shading)
    ax2.pcolormesh(y*lfac,z*lfac,bx0.transpose(),cmap='gist_gray',vmax=3.e3,vmin=-3.e3,shading=shading)
    ax3.pcolormesh(y*lfac,z*lfac,bh0.transpose(),cmap='inferno',vmax=5.e3,vmin=-5.e3,shading=shading)

    ax1.set_title('Emergent intensity')
    ax2.set_title('LoS magnetic field')
    ax3.set_title('Horizontal magnetic field')

    ax1.set_xlabel("Mm")
    ax1.set_ylabel("Mm")

    ax2.set_xlabel("Mm")
    ax3.set_xlabel("Mm")
    
    for ax in [ax2,ax3]:
        ax.tick_params(labelleft=False,left=False)

    bbox_props = dict(boxstyle="round,pad=0.3", fc="white", ec="black", lw=2,alpha=0.9)
    ax1.annotate(text="t="+"{:.2f}".format((t-t0)/3600.)+" [hr]"\
                     ,xy=[0.03,0.03],xycoords="axes fraction"\
                     ,fontsize=18,color='black',bbox=bbox_props)

    if n == n0:
        fig.tight_layout()


    plt.savefig(pngdir+"py"+'{0:08d}'.format(n)+".png")
    plt.pause(0.1)

    if(n != nd_tau):
        plt.clf()
