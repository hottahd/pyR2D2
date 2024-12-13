import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pyR2D2
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
casedir="../figs/"+caseid
os.makedirs(casedir,exist_ok=True)

pngdir="../figs/"+caseid+"/slice/"
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

plt.clf()
plt.close('all')

shading='gouraud'

n0 = 30
nd_tau = n0

fig = plt.figure(1)
for n in range(n0,nd_tau+1):
    print(n)
    #d.read_qq_slice(5,'x',0)
    d.read_qq_tau(n,silent=True)
    #ax = fig.add_subplot(111,projection=ccrs.Orthographic(central_longitude=0,central_latitude=np.pi))
    ax = fig.add_subplot(111,projection='mollweide')
    
    var = 'vx'
    vmax = d.qt_yin[var].max()
    #vmax_yan = d.ql_yan[var].max()
    #vmax = np.array([vmax_yin,vmax_yan]).max()
    vmax = 5.e3
    vmin = -vmax

    ax.pcolormesh(zo_yy[jx_yy//2:jx_yy,:],yo_yy[jx_yy//2:jx_yy,:]-0.5*pi,d.qt_yan[var][jx_yy//2:jx_yy,:]
                  ,vmin=vmin,vmax=vmax,shading=shading,cmap='gray')

    #ax.pcolormesh(zog_yy,yog_yy-0.5*pi,d.ql_yan[var]
    #              ,vmin=vmin,vmax=vmax,shading=shading)
    
    ax.pcolormesh(zo_yy[0:jx_yy//2,:],yo_yy[0:jx_yy//2,:]-0.5*pi,d.qt_yan[var][0:jx_yy//2,:]
                  ,vmin=vmin,vmax=vmax,shading=shading)
    ax.pcolormesh(zz_yy,yy_yy-0.5*pi,d.qt_yin[var],vmin=vmin,vmax=vmax,shading=shading)
    ax.set_xticklabels('')
    ax.set_yticklabels('')

    plt.pause(0.01)
    plt.savefig(pngdir+"py"+'{0:08d}'.format(n)+".png")

    if(n != nd_tau):
        plt.clf()
