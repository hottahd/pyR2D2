import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
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
pngdir="../figs/"+caseid+"/mov_photo_720/"
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

yran = ymax - ymin
xran = min(xmax-xmin,yran)

xsize = 7.2
ysize = 7.2
fig = plt.figure(num=1,figsize=(xsize,ysize))

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

    shading = "auto"

    lfac = 1.e-8
    
    ax1 = fig.add_axes([0,0,1,1])
    
    in0 = d.qt["in"].copy()
    in0s = np.roll(in0,[jx//2-jc,kx//2-kc],axis=[0,1])
    ax1.pcolormesh(y*lfac,z*lfac,in0s.transpose(),cmap='gist_gray',vmax=3.2e10,vmin=1.e10,shading=shading)
    
    plt.pause(0.1)
    plt.savefig(pngdir+"py"+'{0:08d}'.format(n)+".png")

    if(n != nd_tau):
        clf()
