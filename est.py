import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from tqdm import tqdm
import R2D2
import sys

try:
    caseid
except NameError:
    print("input caseid id (3 digit)")
    caseid = 0
    caseid = input()
    caseid = "d"+caseid.zfill(3)

datadir="../run/"+caseid+"/data/"
pngdir="../figs/"+caseid+"/est/"
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

if geometry == 'Cartesian':
    sinyy = 1
    sinyym = 1.
else:
    xx,yy = np.meshgrid(x,y,indexing='ij')
    sinyy = sin(yy)
    sinyym = np.average(sinyy,axis=1)

    xx,yy_flux = np.meshgrid(x_flux,y,indexing='ij')
    sinyy_flux = sin(yy_flux)
    sinyym_flux = np.average(sinyy_flux,axis=1)

#n0 = 4
#nd = n0

#nd = 10

vxrmst = np.zeros((ix,nd-n0+1))
vyrmst = np.zeros((ix,nd-n0+1))
vzrmst = np.zeros((ix,nd-n0+1))

bxmt = np.zeros((ix,jx,nd-n0+1))
bymt = np.zeros((ix,jx,nd-n0+1))
bzmt = np.zeros((ix,jx,nd-n0+1))

bxrmst = np.zeros((ix,nd-n0+1))
byrmst = np.zeros((ix,nd-n0+1))
bzrmst = np.zeros((ix,nd-n0+1))

rormst = np.zeros((ix,nd-n0+1))
sermst = np.zeros((ix,nd-n0+1))
prrmst = np.zeros((ix,nd-n0+1))
termst = np.zeros((ix,nd-n0+1))

romt = np.zeros((ix,nd-n0+1))
semt = np.zeros((ix,nd-n0+1))
prmt = np.zeros((ix,nd-n0+1))
temt = np.zeros((ix,nd-n0+1))

fet = np.zeros((ix+1,nd-n0+1))
fmt = np.zeros((ix+1,nd-n0+1))
fdt = np.zeros((ix+1,nd-n0+1))
fkt = np.zeros((ix+1,nd-n0+1))
frt = np.zeros((ix+1,nd-n0+1))
ftt = np.zeros((ix+1,nd-n0+1))


for n in tqdm(range(n0,nd+1)):
    #print(n)
    ##############################
    # read time
    t = d.read_time(n)
    
    ##############################
    # read value
    d.read_vc(n,silent=True)

    ##############################    
    if geometry == 'Cartesian':
        fstar = lstar/4/np.pi/rstar**2
        fe = np.average(d.vc["fe"],axis=1)/sinyym
        fd = np.average(d.vc["fd"],axis=1)/sinyym
        fk = np.average(d.vc["fk"],axis=1)/sinyym
        fr = np.average(d.vc["fr"],axis=1)/sinyym
        fm = np.average(d.vc["fm"],axis=1)/sinyym

    else:
        fstar = lstar/4/np.pi
        j1 = jx//2 - 512
        j2 = jx//2 + 512
        #fe = np.average(d.vc["fe"]*sinyy_flux,axis=1)/sinyym_flux*x_flux**2
        #fd = np.average(d.vc["fd"]*sinyy_flux,axis=1)/sinyym_flux*x_flux**2

        sinyym_flux0 = np.average(sinyy_flux[:,j1:j2],axis=1)
        
        fe = np.average(d.vc["fe"][:,j1:j2]*sinyy_flux[:,j1:j2],axis=1)/sinyym_flux0*x_flux**2
        fd = np.average(d.vc["fd"][:,j1:j2]*sinyy_flux[:,j1:j2],axis=1)/sinyym_flux0*x_flux**2
        
        fk = np.average(d.vc["fk"]*sinyy_flux,axis=1)/sinyym_flux*x_flux**2
        fm = np.average(d.vc["fm"]*sinyy_flux,axis=1)/sinyym_flux*x_flux**2
        fr = np.average(d.vc["fr"]*sinyy_flux,axis=1)/sinyym_flux#*x_flux**2
    
    xs = rstar - 2.e8
    ds = 2.e7
    sr = 0.5e0*(1.e0 + np.tanh((x_flux-xs)/ds))
    SR, sry = np.meshgrid(sr,y,indexing="ij")
    
    ff = fd*sr + fe*(1.e0-sr)
    #ff = fe
    #ff = fd
    ft = ff + fk + fr + fm

    vxrmst[:,n-n0] = np.sqrt(np.average(d.vc["vxrms"]**2*sinyy,axis=1)/sinyym)
    vyrmst[:,n-n0] = np.sqrt(np.average(d.vc["vyrms"]**2*sinyy,axis=1)/sinyym)
    vzrmst[:,n-n0] = np.sqrt(np.average(d.vc["vzrms"]**2*sinyy,axis=1)/sinyym)

    bxrmst[:,n-n0] = np.sqrt(np.average(d.vc["bxrms"]**2*sinyy,axis=1)/sinyym)
    byrmst[:,n-n0] = np.sqrt(np.average(d.vc["byrms"]**2*sinyy,axis=1)/sinyym)
    bzrmst[:,n-n0] = np.sqrt(np.average(d.vc["bzrms"]**2*sinyy,axis=1)/sinyym)

    bxmt[:,:,n-n0] = d.vc['bxm']
    bymt[:,:,n-n0] = d.vc['bym']
    bzmt[:,:,n-n0] = d.vc['bzm']

    rormst[:,n-n0] = np.sqrt(np.average(d.vc["rorms"]**2*sinyy,axis=1)/sinyym)
    sermst[:,n-n0] = np.sqrt(np.average(d.vc["serms"]**2*sinyy,axis=1)/sinyym)
    prrmst[:,n-n0] = np.sqrt(np.average(d.vc["prrms"]**2*sinyy,axis=1)/sinyym)
    termst[:,n-n0] = np.sqrt(np.average(d.vc["terms"]**2*sinyy,axis=1)/sinyym)

    romt[:,n-n0] = np.average(d.vc["rom"]*sinyy,axis=1)/sinyym
    semt[:,n-n0] = np.average(d.vc["sem"]*sinyy,axis=1)/sinyym
    prmt[:,n-n0] = np.average(d.vc["prm"]*sinyy,axis=1)/sinyym
    temt[:,n-n0] = np.average(d.vc["tem"]*sinyy,axis=1)/sinyym

    fet[:,n-n0] = fe
    fmt[:,n-n0] = fm
    fdt[:,n-n0] = fd
    fkt[:,n-n0] = fk
    frt[:,n-n0] = fr
    ftt[:,n-n0] = ft

    fontsize = 12
    fmax = 2.0
    fmin = -1.0

    fig = plt.figure(num='est',figsize=(12,8))
    ax1 = fig.add_subplot(221)
    ax2 = fig.add_subplot(222)
    ax3 = fig.add_subplot(223)
    ax4 = fig.add_subplot(224)


    #####################
    if geometry == 'Spherical':
        xp = x_flux/rstar
        xlabel = r'$r/R_\odot$'
        xpp = x/rstar
    else:
        xp = (x_flux - rstar)*1.e-8
        xpp = (x - rstar)*1.e-8
        xlabel = r'$x-R_\odot\ \mathrm{[Mm]}$'
        
    ax1.plot(xp,ff/fstar,label=r'$F_\mathrm{e}$',color=R2D2.magenta)
    ax1.plot(xp,fk/fstar,label=r'$F_\mathrm{k}$',color=R2D2.green)
    ax1.plot(xp,fr/fstar,label=r'$F_\mathrm{r}$',color=R2D2.blue)
    ax1.plot(xp,fm/fstar,label=r'$F_\mathrm{r}$',color=R2D2.orange)
    ax1.plot(xp,ft/fstar,label=r'$F_\mathrm{t}$',color=R2D2.ash)

    ax1.hlines(y=1,xmin=xp.min(),xmax=xp.max(),linestyle='--',color=R2D2.ash)
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel("$F/F_{\odot}$")
    ax1.set_title("Energy fluxes")
    #ax1.legend()

    #####################
    vxrms = vxrmst[:,n-n0]
    vhrms = np.sqrt(vyrmst[:,n-n0]**2 + vzrmst[:,n-n0]**2)
    vrms = np.sqrt(vxrms**2 + vhrms**2)
    ax2.plot(xpp,vxrms*1.e-5,label=r'$v_{x\mathrm{(rms)}}$',color=R2D2.blue)
    ax2.plot(xpp,vhrms*1.e-5,label=r'$v_\mathrm{h(rms)}$',color=R2D2.magenta)
    ax2.plot(xpp,vrms*1.e-5,label=r'$v_\mathrm{(rms)}$',color=R2D2.green)
    ax2.set_xlabel(xlabel)
    ax2.set_ylabel(r"velocities [km/s]")
    ax2.set_label('RMS velocities')
    if deep_flag == 0:
        ax2.set_yscale('log')
    ax2.legend()

    #####################
    bxrms = bxrmst[:,n-n0]
    bhrms = np.sqrt(byrmst[:,n-n0]**2 + bzrmst[:,n-n0]**2)
    brms = np.sqrt(bxrms**2 + bhrms**2)
    ax3.plot(xpp,bxrms,label=r'$B_{x\mathrm{(rms)}}$',color=R2D2.blue)
    ax3.plot(xpp,bhrms,label=r'$B_\mathrm{h(rms)}$',color=R2D2.magenta)
    ax3.plot(xpp,brms,label=r'$B_\mathrm{(rms)}$',color=R2D2.green)
    ax3.set_xlabel(xlabel)
    ax3.set_ylabel(r"Magnetic field [G]")
    ax3.set_label('RMS magnetic field')
    ax3.legend()

    #####################
    ax4.plot(x/rstar,semt[:,n-n0]+se0,color=R2D2.blue)
        
    if n == n0:
        fig.tight_layout()

    ax3.annotate(text="t="+"{:.2f}".format(t/3600./24.)+" [day]"\
                     ,xy=[0.01,0.01],xycoords="figure fraction",fontsize=18)

    plt.pause(0.001)
    plt.savefig(pngdir+"py"+'{0:08d}'.format(n)+".png")

    
    if n != nd:
        plt.clf() # clear figure

    # loop end
    ###############################################################################
    ###############################################################################
    ###############################################################################

vxrms = np.sqrt(np.average(vxrmst**2,axis=1))
vyrms = np.sqrt(np.average(vyrmst**2,axis=1))
vzrms = np.sqrt(np.average(vzrmst**2,axis=1))

bxrms = np.sqrt(np.average(bxrmst**2,axis=1))
byrms = np.sqrt(np.average(byrmst**2,axis=1))
bzrms = np.sqrt(np.average(bzrmst**2,axis=1))

rorms = np.sqrt(np.average(rormst**2,axis=1))
serms = np.sqrt(np.average(sermst**2,axis=1))
prrms = np.sqrt(np.average(prrmst**2,axis=1))
terms = np.sqrt(np.average(termst**2,axis=1))

rom = np.average(romt,axis=1)
sem = np.average(semt,axis=1)
prm = np.average(prmt,axis=1)
tem = np.average(temt,axis=1)

fe = np.average(fet,axis=1)
fd = np.average(fdt,axis=1)

ff = fd*sr + fe*(1.e0-sr)
#ff = fd

fk = np.average(fkt,axis=1)
fr = np.average(frt,axis=1)
ft = ff + fk + fr + fm

np.savez(d.p['datadir']+"est.npz"\
             ,x=x,y=y,z=z,rstar=rstar\
             ,ro0=ro0,pr0=pr0,te0=te0,se0=se0\
             ,vxrms=vxrms,vyrms=vyrms,vzrms=vzrms\
             ,bxrms=bxrms,byrms=byrms,bzrms=bzrms\
             ,rorms=rorms,prrms=prrms,serms=serms,terms=terms\
             ,rom=rom,prm=prm,sem=sem,tem=tem\
             ,ff=ff,fk=fk,fr=fr,ft=ft\
             )
         
fig2 = plt.figure(num=100,figsize=(12,5))
ax23 = fig2.add_subplot(121)
ax24 = fig2.add_subplot(122)
ax23.plot(xp,ff/fstar,color=R2D2.magenta,label="$F_\mathrm{e}$")
ax23.plot(xp,fk/fstar,color=R2D2.green,label="$F_\mathrm{k}$")
ax23.plot(xp,fr/fstar,color=R2D2.blue,label="$F_\mathrm{r}$")
ax23.plot(xp,ft/fstar,color=R2D2.orange,label="$F_\mathrm{t}$")
ax23.plot(xp,fm/fstar,color=R2D2.ash,label="$F_\mathrm{m}$")
#ax23.set_xlim(xmin/rsun,xmax/rsun)
ax23.set_ylim(fmin,fmax)
ax23.set_xlabel(xlabel)
ax23.set_ylabel("$F/F_{\odot}$")
ax23.set_title("Full convection zone")
ax23.legend(loc='upper left',prop={'size': 15})
ax23.annotate(text="t="+"{:.2f}".format(t/3600./24.)+" [day]"\
                 ,xy=[0.01,0.01],xycoords="figure fraction",fontsize=18)

ax23.hlines(y=1,xmin=xp.min(),xmax=xp.max(),linestyle='--',color='black')

x_flux_c = (x_flux - rstar)*1.e-8

ax24.plot(x_flux_c,ff/fstar,color=R2D2.magenta)
ax24.plot(x_flux_c,fk/fstar,color=R2D2.green)
ax24.plot(x_flux_c,fr/fstar,color=R2D2.blue)
ax24.plot(x_flux_c,ft/fstar,color=R2D2.orange)
ax24.plot(x_flux_c,fm/fstar,color=R2D2.ash,label="$F_\mathrm{m}$")
ax24.hlines(y=1,xmin=x_flux_c.min(),xmax=x_flux_c.max(),linestyle='--',color='black')
ax24.set_xlim(-10,1)
ax24.set_ylim(fmin,fmax)
ax24.set_xlabel("$x - R_{\odot} \ [\mathrm{Mm}]$")
ax24.set_ylabel("$F/F_{\odot}$")
ax24.set_title("Around photosphere")
fig2.tight_layout()
plt.pause(0.001)
plt.ion()
    
