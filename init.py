import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
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
casedir="../figs/"+caseid
os.makedirs(casedir,exist_ok=True)

d = R2D2.R2D2_data(datadir)
for key in d.p:
    exec('%s = %s%s%s' % (key, 'd.p["',key,'"]'))

try:
    n0
except NameError:
    n0 = 0
if  n0 > d.p["nd"]:
    n0 = d.p["nd"]

print('### calculation domain ###')
print('xmax - rsun = ', '{:6.2f}'.format((xmax - rsun)*1.e-8),'[Mm], xmin - rsun = ', '{:.2f}'.format((xmin - rsun)*1.e-8),'[Mm]')
print('ymax        = ', '{:6.2f}'.format(ymax*1.e-8)       ,'[Mm], ymin        = ', '{:.2f}'.format(ymin*1.e-8),'[Mm]' )
print('zmax        = ', '{:6.2f}'.format(zmax*1.e-8)       ,'[Mm], zmin        = ', '{:.2f}'.format(zmin*1.e-8),'[Mm]' )

print('')
print('### number of grid ###')
print('(ix,jx,kx)=(',ix,',',jx,',',kx,')')

print('')
print('### calculation time ###')
print('time step (nd) =',nd)
t = d.read_time(nd)
print('time =','{:.2f}'.format(t/3600),' [hour]')
