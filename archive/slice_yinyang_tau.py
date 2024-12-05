import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import R2D2
import sys
import os
import cartopy.crs as ccrs
import shapely.geometry as sgeom

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

shading='flat'

n0 = 80
nd_tau = 80

n0 = 40
nd_tau = 40

fig = plt.figure(100,figsize=(18,9))
for n in range(n0,nd_tau+1):
    print(n)
    #d.read_qq_slice(5,'x',0)
    #d.read_qq_slice(3,'x',n,silent=True)
    d.read_qq_tau(n,silent=True)    

    var = 'in'
    qq_yin = d.qt_yin[var]
    qq_yan = d.qt_yan[var]

    #qq_yin = d.ql_yin[var][margin:jxg_yy-margin,margin:kxg_yy-margin]
    #qq_yan = d.ql_yan[var][margin:jxg_yy-margin,margin:kxg_yy-margin]
    
    qqm = np.array([qq_yin,qq_yan]).mean()

    qq_yin = qq_yin - qqm
    qq_yan = qq_yan - qqm
    
    ax1 = fig.add_subplot(121,projection=ccrs.Orthographic(central_longitude=0.0,central_latitude=10.0))
    ax2 = fig.add_subplot(122,aspect='equal')
    #ax1 = fig.add_subplot(111,projection=ccrs.Mollweide())
    
    #vmax = d.ql_yin[var].max()
    #vmax_yan = d.ql_yan[var].max()
    
    vmax = np.array([qq_yin,qq_yan]).max()*0.5
    vmin = -vmax
    #vmax = 5.e3
    #vmin = -vmax

    rad2deg = 180/pi
    
    ax1.pcolormesh((zo_yy[jx_yy//2:jx_yy,:])*rad2deg,(yo_yy[jx_yy//2:jx_yy,:]-0.5*pi)*rad2deg,qq_yan[jx_yy//2:jx_yy,:]
                ,transform=ccrs.PlateCarree(),vmin=vmin,vmax=vmax,shading=shading,cmap='gray')
    ax1.pcolormesh(zo_yy[0:jx_yy//2,:]*rad2deg,(yo_yy[0:jx_yy//2,:]-0.5*pi)*rad2deg,qq_yan[0:jx_yy//2,:] \
                  ,transform=ccrs.PlateCarree(),vmin=vmin,vmax=vmax,shading=shading)
    ax1.pcolormesh((zz_yy)*rad2deg,(yy_yy-0.5*pi)*rad2deg,qq_yin,transform=ccrs.PlateCarree(),vmin=vmin,vmax=vmax,shading=shading)

    #box = sgeom.box(minx=30,maxx=50,miny=40,maxy=60)
    #ax1.add_geometries([box],ccrs.PlateCarree(),facecolor='none',edgecolor='black')

    yt_yy = (y_yy - 0.5*pi + 0.5*(y_yy[1] - y_yy[0]) )*rad2deg
    zt_yy = (z_yy + 0.5*(z_yy[1] - z_yy[0]))*rad2deg

    y0 = -5
    y1 = +5
    z0 = -5
    z1 = +5
    
    j0 = np.argmin(np.abs(yt_yy - y0))
    j1 = np.argmin(np.abs(yt_yy - y1))

    k0 = np.argmin(np.abs(zt_yy - z0))
    k1 = np.argmin(np.abs(zt_yy - z1))
    
    #j0 = jx_yy//2 - jx_yy//10
    #j1 = jx_yy//2 + jx_yy//10
    #k0 = kx_yy//2 - kx_yy//10
    #k1 = kx_yy//2 + kx_yy//10
    
    ax2.pcolormesh((zz_yy[j0:j1,k0:k1])*rad2deg,(yy_yy[j0:j1,k0:k1]-0.5*pi)*rad2deg,qq_yin[j0:j1,k0:k1],vmin=vmin,vmax=vmax,shading=shading)
    
    #ax.pcolormesh(zog_yy[jxg_yy//2:jxg_yy,:],yog_yy[jxg_yy//2:jxg_yy,:]-0.5*pi,d.ql_yan[var][jxg_yy//2:jxg_yy,:]
    #              ,vmin=vmin,vmax=vmax,shading=shading,cmap='gray')

    
    #ax.pcolormesh(zog_yy,yog_yy-0.5*pi,d.ql_yan[var]
    #              ,vmin=vmin,vmax=vmax,shading=shading)

    ax1.plot(zt_yy[k0:k1],y0*np.ones(k1-k0),transform=ccrs.PlateCarree(),color='black')
    ax1.plot(zt_yy[k0:k1],y1*np.ones(k1-k0),transform=ccrs.PlateCarree(),color='black')
    ax1.plot(z0*np.ones(j1-j0),yt_yy[j0:j1],transform=ccrs.PlateCarree(),color='black')
    ax1.plot(z1*np.ones(j1-j0),yt_yy[j0:j1],transform=ccrs.PlateCarree(),color='black')

    #if n == n0:
    #    plt.tight_layout()
    plt.subplots_adjust(left=0.02,right=0.98,bottom=0.02,top=0.98,wspace=0.1,hspace=0.1)
    
    #ax.pcolormesh(zog_yy[0:jxg_yy//2,:],yog_yy[0:jxg_yy//2,:]-0.5*pi,d.ql_yan[var][0:jxg_yy//2,:]              
    #              ,vmin=vmin,vmax=vmax,shading=shading)
    #ax.set_xticklabels('')
    #ax.set_yticklabels('')

    plt.pause(0.01)
    plt.savefig(pngdir+"py"+'{0:08d}'.format(n)+".png")

    if(n != nd_tau):
        plt.clf()
