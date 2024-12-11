def add_color_bar_2x2(fig,axes,ims,ysize):
    '''
    put color bar
    fig: figure instance
    axes: 
    ims: output of pcolormesh
    '''
    # add color bar
    box3 = axes[2].get_position().bounds
    wid3, hei3 = box3[2], box3[3]
    
    color_bar_width = 0.12
    color_bar_height = 0.08/ysize
    
    for ax, im in zip(axes,ims):
        box = ax.get_position().bounds
        cax = fig.add_axes([box[0]+ box[2] - color_bar_width
                        ,box[1] + box[3] + color_bar_height*4
                        ,color_bar_width
                        ,color_bar_height
                        ])
        cbar = fig.colorbar(im,cax=cax,orientation='horizontal')
        cbar.ax.tick_params(labelsize=12)

def mov_cartesian_photo_2x2(data,t,vls,tu_height,vmaxs,vmins,titles,cmaps=['inferno','gray','inferno','gray'],tight_layout_flag=True):
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec

    zran = data.zmax - data.zmin
    yran = data.ymax - data.ymin
    xran = min(data.xmax - data.xmin, zran)

    xsize = 18
    ysize = xsize*(yran + xran)/2/zran
    ysize = ysize + min(xsize/ysize/1.5,2)
    ysize_limit = 11
    if ysize > ysize_limit:
        xsize = xsize/ysize*ysize_limit
        ysize = ysize_limit
    fig = plt.figure(num=1,figsize=(xsize,ysize))

    grid = GridSpec(2,2,height_ratios=[yran,xran])    
    
    ax1 = fig.add_subplot(grid[0,0],aspect='equal')
    ax2 = fig.add_subplot(grid[0,1],aspect='equal')
    ax3 = fig.add_subplot(grid[1,0],aspect='equal')
    ax4 = fig.add_subplot(grid[1,1],aspect='equal')
        
    shading = "auto"
    lfac = 1.e-8 # unit is Mm
    
    ims = []
    axes = [ax1,ax2,ax3,ax4]
    # pcolormesh
    for ax, vl, cmap, vmax, vmin in zip(axes[:2],vls[:2],cmaps[:2],vmaxs[:2],vmins[:2]):
        ims.append(ax.pcolormesh(data.z*lfac, data.y*lfac, vl, cmap=cmap, \
                vmax=vmax, vmin=vmin, shading=shading))
    
    for ax, vl, cmap, vmax, vmin in zip(axes[2:],vls[2:],cmaps[2:],vmaxs[2:],vmins[2:]):
        ims.append(ax.pcolormesh(data.z*lfac,(data.x - data.rstar)*lfac, vl \
                ,cmap=cmap, vmin=vmin, vmax=vmax, shading=shading))
        # to deal with older version of R2D2
        if tu_height.max() > 0.8*data.rstar:
            tu_height = tu_height - data.rstar
        ax.plot(data.z*lfac, tu_height*lfac,color="w")
    
    for ax in [ax1,ax2]:
        ax.tick_params(labelbottom=False)
    for ax in [ax2,ax4]:
        ax.tick_params(labelleft=False)
    for ax in [ax3,ax4]:
        ax.set_xlabel('$z$ [Mm]')
    ax1.set_ylabel('$y$ [Mm]')
    ax3.set_ylabel('$x$ [Mm]')
    
    for ax in [ax3,ax4]:
        ax.set_ylim((max(data.xmax-yran,data.xmin) - data.rstar)*lfac,
                    (data.xmax - data.rstar)*lfac)
    
    for ax, title in zip([ax1,ax2,ax3,ax4],titles):
        title = ax.set_title(title)
        title.set_position([0.01,1.02])
        title.set_ha('left')
        
    if tight_layout_flag:
        fig.tight_layout()
        
    ax3.annotate(text="$t="+"{:.2f}".format((t)/60/60)+"~\mathrm{[hour]}$"\
                     ,xy=[0.01,0.02],xycoords="figure fraction"\
                     ,color='black',fontsize=20)
    
    # add color bar
    add_color_bar_2x2(fig,axes,ims,ysize)

def mov_spherical_2x2(data,t,vls,vmaxs,vmins,titles,cmaps=['inferno','gray','inferno','gray'],tight_layout_flag=True):
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec
        
    xsize = 16
    ysize = 9
    fig = plt.figure(num=1,figsize=(xsize,ysize))
        
    ax1 = fig.add_subplot(221,projection='mollweide')
    ax2 = fig.add_subplot(222,projection='mollweide')
    ax3 = fig.add_subplot(223,aspect='equal')
    ax4 = fig.add_subplot(224,aspect='equal')
    
    shading = "auto"    
    ims = []
    axes = [ax1,ax2,ax3,ax4]
    
    lfac = 1/data.rstar
    
    for ax, title in zip(axes,titles):
        title = ax.set_title(title)
        title.set_position([0.01,1.02])
        title.set_ha('left')
        
                
    for ax, vl, cmap, vmax, vmin in zip(axes[:2],vls[:2],cmaps[:2],vmaxs[:2],vmins[:2]):
        ims.append(ax.pcolormesh(data.zz, data.yy ,vl \
                ,shading=shading, cmap=cmap, vmax=vmax, vmin=vmin))
        ax.set_xticklabels('')
        ax.set_yticklabels('')
    
    for ax, vl, cmap, vmax, vmin in zip(axes[2:],vls[2:],cmaps[2:],vmaxs[2:],vmins[2:]):
        ims.append(ax.pcolormesh(data.XX.T*lfac, data.YY.T*lfac, vl.T
                                 ,shading=shading, cmap=cmap, vmax=vmax, vmin=vmin))
        ax.set_xlabel(r'$z/R_*$')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        
    ax3.set_ylabel(r'$x/R_*$')
        
    if tight_layout_flag:
        fig.tight_layout()
        
    ax3.annotate(text="$t="+"{:.2f}".format((t)/60/60/24)+"~\mathrm{[day]}$"\
                    ,xy=[0.01,0.02],xycoords="figure fraction"\
                    ,color='black',fontsize=20)
    
    # add color bar
    add_color_bar_2x2(fig,axes,ims,ysize)
    
def mov_yinyang_2(data,t,vls,vmaxs,vmins,titles,
                  cmaps=['inferno','gray'],
                  central_longitude=0,
                  central_latitude=30,
                  tight_layout_flag=True):
    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs
    import numpy as np

    data.YinYangSet()
    
    xsize = 16
    ysize = 9
    fig = plt.figure(num=1,figsize=(xsize,ysize))
            
    ax1 = fig.add_subplot(121,
                          projection=ccrs.Orthographic(central_longitude=central_longitude,central_latitude=central_latitude))
    ax2 = fig.add_subplot(122,
                          projection=ccrs.Orthographic(central_longitude=central_longitude,central_latitude=central_latitude))
    
    rad2deg = 180/np.pi
    axes = [ax1,ax2]
    ims = []
    for vl, ax, cmap, vmax, vmin, title in zip(vls, axes, cmaps, vmaxs, vmins, titles):
        for z, y, YinYang in zip(['Zog_yy', 'Zg_yy'],['Yog_yy', 'Yg_yy'],['Yan', 'Yin']):
            ax.pcolormesh(data.p.__dict__[z]*rad2deg, \
                (data.p.__dict__[y]-0.5*np.pi)*rad2deg, vl[YinYang], \
                transform=ccrs.PlateCarree(), cmap=cmap, vmax=vmax, vmin=vmin)
        ims.append(ax.collections[0])
        ax.set_title(title)
        ax.set_xticklabels('')
        ax.set_yticklabels('')

    if tight_layout_flag:
        fig.tight_layout()
        
    # add color bar
    
    color_bar_width = 0.12
    color_bar_height = 0.08/ysize
    
    for ax, im in zip(axes,ims):
        box = ax.get_position().bounds
        cax = fig.add_axes([box[0]+ box[2] - color_bar_width
                        ,box[1] + color_bar_height
                        ,color_bar_width
                        ,color_bar_height
                        ])
        cbar = fig.colorbar(im,cax=cax,orientation='horizontal')
        cbar.ax.tick_params(labelsize=12)
        
    ax1.annotate(text="$t="+"{:.2f}".format((t)/60/60/24)+"~\mathrm{[day]}$"\
            ,xy=[0.01,0.02],xycoords="figure fraction"\
            ,color='black',fontsize=25)