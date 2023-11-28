import numpy as np
import matplotlib.pyplot as plt
import R2D2
import sys
import os
import cartopy.crs as ccrs

try:
    caseid
except NameError:
    print("input caseid id (3 digit)")
    caseid = 0
    caseid = input()
    caseid = "d"+caseid.zfill(3)
    
datadir="../run/"+caseid+"/data/"
pngdir="../figs/"+caseid+"/png/"
os.makedirs(pngdir,exist_ok=True)

d = R2D2.R2D2_data(datadir)
for key in d.p:
    exec('%s = %s%s%s' % (key, 'd.p["',key,'"]'))
    
try:
    n0
except NameError:
    n0 = 0
if  n0 > d.p["nd"]:
    n0 = d.p["nd"]

print("Maximum time step= ",nd," time ="\
          ,dtout*float(nd)/3600./24.," [day]")

plt.close('all')

X, Y = np.meshgrid(x,y,indexing='ij')
SINY = np.sin(Y)
SINYM = np.sin(Y).sum(axis=1)

plt.rcParams['font.size'] = 16
# read time
t0 = d.read_time(0,silent=True)

yran = ymax - ymin
xran = min(xmax-xmin,yran)

xsize = 20
ysize = 12
fig = plt.figure(num=1,figsize=(xsize,ysize))

#grid = GridSpec(2,2,height_ratios=[yran,xran])

xe = np.zeros(ix+1)
ye = np.zeros(jx+1)

xe[0] = xmin
xe[ix] = xmax
ye[0] = 0.e0
ye[jx] = np.pi
for i in range(1,ix):
    xe[i] = xe[i-1] + 2*(x[i-1] - xe[i-1])

for j in range(1,jx):
    ye[j] = ye[j-1] + 2*(y[j-1] - ye[j-1])

RA, TH = np.meshgrid(x,y,indexing='ij')

RAE, THE = np.meshgrid(xe,ye,indexing='ij')
XX, YY = RA*np.cos(THE), RAE*np.sin(THE)

zz,yy = np.meshgrid(z,y-0.5*np.pi)

yc = y/np.pi*180 - 90
zc = z/np.pi*180

YC, ZC = np.meshgrid(yc,zc,indexing='ij')

for n in range(n0,nd+1):
#for n in range(0,1):
    print(n)
    ##############################
    # read time
    t = d.read_time(n,silent=True)
    
    ##############################
    # read value

    d.read_vc(n,silent=True)
    ##############################

    shading = "flat"
    #shading = "groroud"

    ax1 = fig.add_subplot(221,projection='mollweide')
    ax2 = fig.add_subplot(222,projection='mollweide')
    #ax1 = fig.add_subplot(221,aspect='equal')
    #ax2 = fig.add_subplot(222,aspect='equal')
    ax3 = fig.add_subplot(223,aspect='equal')
    ax4 = fig.add_subplot(224,aspect='equal')

    if xmax > rstar:
        d.read_qq_tau(n*int(ifac),silent=True)
        ax1.pcolormesh(zz,yy,d.qt['vx'],shading='auto')
        ax2.pcolormesh(zz,yy,d.qt['bx'],shading='auto')
    else:
        d.read_qq_select(xmax,n,silent=True)
        vx = d.qs['vx']
        bx = d.qs['bx']
        vxrms = np.sqrt((vx**2).mean())
        vmax = 2*vxrms
        vmin = -vmax
        ax1.pcolormesh(zz,yy,vx,shading='auto',cmap='inferno',vmax=vmax,vmin=vmin)
        bxrms = np.sqrt((bx**2).mean())
        vmax = 2*bxrms
        vmin = -vmax
        ax2.pcolormesh(zz,yy,bx,shading='auto',cmap='gray',vmax=vmax,vmin=vmin)
        for ax in [ax1,ax2]:
            ax.set_xticklabels('')
            ax.set_yticklabels('')
        
    sem, tmp   = np.meshgrid((d.vc['sem']*SINY).sum(axis=1)/SINYM,y,indexing='ij')
    serms, tmp = np.meshgrid(np.sqrt((d.vc['serms']**2*SINY).sum(axis=1)/SINYM),y,indexing='ij')
    bbp = np.sqrt(d.vc['bx_xy']**2 + d.vc['by_xy']**2 + d.vc['bz_xy']**2)
    om = d.vc['vzm']/RA/np.sin(TH)
    bzm = d.vc['bzm']
    
    if serms.max() != 0:
        se_plot = (d.vc['se_xy']-sem)/serms
    else:
        se_plot = np.zeros((ix,jx))

    lfac = 1/rstar
    ax3.pcolormesh(XX.T*lfac,YY.T*lfac,se_plot.T,vmin=-2.,vmax=2.,shading='auto')
    ax4.pcolormesh(XX.T*lfac,YY.T*lfac,bzm.T,shading='auto',vmin=-8000.,vmax=8000)
    
    if(n == n0):
        fig.tight_layout()
    
    plt.pause(0.1)
    plt.savefig(pngdir+"py"+'{0:08d}'.format(n)+".png")

    if(n != nd):
        plt.clf()
