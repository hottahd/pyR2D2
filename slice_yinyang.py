import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
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
casedir="../figs/"+caseid
os.makedirs(casedir,exist_ok=True)

pngdir="../figs/"+caseid+"/slice/"
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

plt.clf()
plt.close('all')

shading='gouraud'

n0 = 30
nd_tau = n0

fig = plt.figure(1)
for n in range(n0,nd_tau+1):
    print(n)
    #d.read_qq_slice(5,'x',0)
    d.read_qq_slice(5,'x',n,silent=True)
    ax = fig.add_subplot(111,projection=ccrs.Orthographic(central_longitude=0,central_latitude=np.pi))
    
    var = 'vx'
    vmax = d.ql_yin[var].max()
    #vmax_yan = d.ql_yan[var].max()
    #vmax = np.array([vmax_yin,vmax_yan]).max()
    vmax = 5.e3
    vmin = -vmax

    ax.pcolormesh(zog_yy[jxg_yy//2:jxg_yy,:],yog_yy[jxg_yy//2:jxg_yy,:]-0.5*pi,d.ql_yan[var][jxg_yy//2:jxg_yy,:]
                  ,vmin=vmin,vmax=vmax,shading=shading,cmap='gray')

    #ax.pcolormesh(zog_yy,yog_yy-0.5*pi,d.ql_yan[var]
    #              ,vmin=vmin,vmax=vmax,shading=shading)
    
    ax.pcolormesh(zog_yy[0:jxg_yy//2,:],yog_yy[0:jxg_yy//2,:]-0.5*pi,d.ql_yan[var][0:jxg_yy//2,:]              
                  ,vmin=vmin,vmax=vmax,shading=shading)
    ax.pcolormesh(zzg_yy,yyg_yy-0.5*pi,d.ql_yin[var],vmin=vmin,vmax=vmax,shading=shading)
    ax.set_xticklabels('')
    ax.set_yticklabels('')

    plt.pause(0.01)
    plt.savefig(pngdir+"py"+'{0:08d}'.format(n)+".png")

    if(n != nd_tau):
        plt.clf()
