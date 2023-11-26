import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import R2D2
import sys,os
    
pngdir="../figs/comp/mov_star_comp_sunspot/"
os.makedirs(pngdir,exist_ok=True)

d02 = R2D2.R2D2_data('../run/d007/data/')
d06 = R2D2.R2D2_data('../run/d008/data/')
d10 = R2D2.R2D2_data('../run/d009/data/')
    
plt.close('all')
plt.rcParams['font.size'] = 22

xsize = 24
ysize = 8
fig = plt.figure('mov_photo',figsize=(xsize,ysize))

n0 = 0
n1 = 4001
for n in tqdm(range(n0,n1)):
    #print(n)
    ##############################
    # read time
    t = d02.read_time(n,tau=True,silent=True)
        
    ##############################
    # read time
    d02.read_qq_tau(n,silent=True)
    d06.read_qq_tau(n,silent=True)
    d10.read_qq_tau(n,silent=True)

    ##############################

    shading = "auto"

    lf = 1.e-8    
    ax1 = fig.add_subplot(131,aspect='equal')
    ax2 = fig.add_subplot(132,aspect='equal')
    ax3 = fig.add_subplot(133,aspect='equal')

    in02_0 = 2.218e9
    in06_0 = 5.83e9
    in10_0 = 2.40e10

    im1 = ax1.pcolormesh(d02.p['y']*lf,d02.p['z']*lf,d02.qt['in'].T/in02_0,vmin=0.65,vmax=1.05)
    im2 = ax2.pcolormesh(d06.p['y']*lf,d06.p['z']*lf,d06.qt['in'].T/in06_0,vmin=0.32,vmax=1.15)
    im3 = ax3.pcolormesh(d10.p['y']*lf,d10.p['z']*lf,d10.qt['in'].T/in10_0,vmin=0.1,vmax=1.4)

    ax1.set_ylabel('Mm')
    titles = [
        r'$0.2M_\odot:~I/I_0=0.65-1.05$',
        r'$0.6M_\odot:~I/I_0=0.32-1.15$',
        r'$M_\odot:~I/I_0=0.1-1.4$'
    ]
    for ax,title in zip([ax1,ax2,ax3],titles):
        ax.set_xlabel('Mm')
        ax.set_title(title)

    if(n == n0):
        fig.tight_layout()

        
    ax1.annotate(text="$t$="+"{:.2f}".format(t/60/60)+" [hour]"\
                     ,xy=[0.01,0.01],xycoords="figure fraction",fontsize=22 \
                     ,color='black')#,path_effects=[withStroke(foreground='black',linewidth=3)])
        
    plt.pause(0.1)
    plt.savefig(pngdir+"py"+'{0:08d}'.format(n-n0)+".png")

    if(n != n1 - 1):
        plt.clf()
