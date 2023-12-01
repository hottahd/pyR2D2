import R2D2
import os

n = 60
ixc = 1024
jxc = 1024
kxc = 1024
var = 'se'

READ    = True
CONVERT = True
WRITE   = True

file = 'd011.vtk'

print('### Time step =',n)
if READ:
    print('- Read data')
    d.read_qq(n,var,silent=True)
    d.read_vc(n,silent=True)

    qqm0  ,tmp = np.meshgrid(d.vc[var+'m']  ,z,indexing='ij')
    qqrms0,tmp = np.meshgrid(d.vc[var+'rms'],z,indexing='ij')

    qqm   = qqm0.reshape((ix,jx,kx))
    qqrms = qqrms0.reshape((ix,jx,kx))

    qqs = (d.qq[var] - qqm)/qqrms

if CONVERT:
    print('- Convert Spherical to Cartesian ###')
    qqc,xc,yc,zc = R2D2.spherical2cartesian(x/rstar,y,z,qqs,ixc,jxc,kxc)

if WRITE:
    print('- Write data in VTK format')
    R2D2.write_3D(qqc,xc,yc,zc,file,var)